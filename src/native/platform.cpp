#include "platform.h"

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <limits.h>
#include <memory>
#include <signal.h>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <bzlib.h>
#include <zlib.h>

namespace native {
namespace {
volatile sig_atomic_t active_child = -1;
volatile sig_atomic_t interrupted = 0;
void signal_handler(int signal) {
  interrupted = signal;
  if (active_child > 0) kill(-active_child, SIGTERM);
}
std::string failure(const std::string &operation, const std::string &path) {
  return operation + " " + path + ": " + std::strerror(errno);
}
int wait_child(pid_t pid) {
  int status;
  while (waitpid(pid, &status, 0) < 0) {
    if (errno != EINTR) throw Error("Cannot wait for child process");
  }
  return WIFEXITED(status) ? WEXITSTATUS(status) : 128 + WTERMSIG(status);
}
void write_all(int fd, const char *data, size_t length) {
  while (length) {
    ssize_t count = write(fd, data, length);
    if (count < 0 && errno == EINTR) continue;
    if (count <= 0) throw Error("Cannot write compressed input/output: " + std::string(std::strerror(errno)));
    data += count;
    length -= count;
  }
}
void fastq_records(const std::string &path, const std::function<void(const std::string &)> &emit) {
  std::ifstream input(path, std::ios::binary);
  if (!input) throw Error("Cannot open FLASH output: " + path);
  std::string header, sequence, separator, quality, buffer;
  buffer.reserve(1024 * 1024 + 65536);
  while (std::getline(input, header)) {
    if (interrupted) throw Error("Interrupted");
    if (!std::getline(input, sequence) || !std::getline(input, separator) || !std::getline(input, quality))
      throw Error("Truncated merged FASTQ: " + path);
    for (auto *line : {&header, &sequence, &separator, &quality})
      if (!line->empty() && line->back() == '\r') line->pop_back();
    if (header.empty() || header[0] != '@' || separator.empty() || separator[0] != '+' ||
        sequence.empty() || sequence.size() != quality.size())
      throw Error("Invalid merged FASTQ: " + path);
    buffer += '>';
    buffer.append(header, 1, std::string::npos);
    buffer += '\n'; buffer += sequence; buffer += '\n';
    if (buffer.size() >= 1024 * 1024) { emit(buffer); buffer.clear(); }
  }
  if (input.bad()) throw Error("Cannot read merged FASTQ: " + path);
  if (!buffer.empty()) emit(buffer);
}
}  // namespace

bool exists(const std::string &path) { struct stat info; return stat(path.c_str(), &info) == 0; }
bool directory(const std::string &path) {
  struct stat info; return stat(path.c_str(), &info) == 0 && S_ISDIR(info.st_mode);
}
uint64_t size(const std::string &path) {
  struct stat info;
  if (stat(path.c_str(), &info)) throw Error(failure("Cannot stat", path));
  return info.st_size;
}
std::string parent(const std::string &path) {
  auto position = path.find_last_of('/');
  return position == std::string::npos ? "." : position == 0 ? "/" : path.substr(0, position);
}
std::string base(const std::string &path) { return path.substr(path.find_last_of('/') + 1); }
bool ends_with(const std::string &value, const std::string &suffix) {
  return value.size() >= suffix.size() && value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}
std::string absolute(const std::string &input) {
  if (input.empty()) throw Error("Empty path");
  std::string path = input;
  if (path[0] == '~' && (path.size() == 1 || path[1] == '/')) {
    const char *home = std::getenv("HOME");
    if (!home) throw Error("HOME is not set");
    path = std::string(home) + path.substr(1);
  }
  if (path[0] != '/') {
    char cwd[PATH_MAX];
    if (!getcwd(cwd, sizeof(cwd))) throw Error("Cannot determine working directory");
    path = std::string(cwd) + "/" + path;
  }
  std::vector<std::string> parts;
  std::stringstream stream(path); std::string part;
  while (std::getline(stream, part, '/')) {
    if (part.empty() || part == ".") continue;
    if (part == "..") { if (!parts.empty()) parts.pop_back(); }
    else parts.push_back(part);
  }
  std::string result;
  for (const auto &piece : parts) result += "/" + piece;
  return result.empty() ? "/" : result;
}
void mkdirs(const std::string &path) {
  if (directory(path)) return;
  if (path.empty() || path == "/" || path == ".") throw Error("Cannot create directory: " + path);
  mkdirs(parent(path));
  if (mkdir(path.c_str(), 0777) && errno != EEXIST) throw Error(failure("Cannot create directory", path));
  if (!directory(path)) throw Error("Not a directory: " + path);
}
void remove_file(const std::string &path) {
  if (unlink(path.c_str()) && errno != ENOENT) throw Error(failure("Cannot remove", path));
}
std::vector<std::string> entries(const std::string &path) {
  std::vector<std::string> result;
  DIR *handle = opendir(path.c_str());
  if (!handle) throw Error(failure("Cannot list", path));
  while (auto *entry = readdir(handle)) {
    std::string name = entry->d_name;
    if (name != "." && name != "..") result.push_back(name);
  }
  closedir(handle);
  std::sort(result.begin(), result.end());
  return result;
}
void remove_tree(const std::string &path) {
  if (path.empty() || absolute(path) == "/") throw Error("Refuse to remove root directory");
  struct stat info;
  if (lstat(path.c_str(), &info)) {
    if (errno == ENOENT) return;
    throw Error(failure("Cannot inspect", path));
  }
  if (!S_ISDIR(info.st_mode)) { remove_file(path); return; }
  for (const auto &name : entries(path)) remove_tree(path + "/" + name);
  if (rmdir(path.c_str())) throw Error(failure("Cannot remove directory", path));
}
void replace(const std::string &source, const std::string &dest) {
  if (rename(source.c_str(), dest.c_str())) throw Error(failure("Cannot publish", dest));
}
void copy(const std::string &source, const std::string &dest) {
  std::ifstream input(source, std::ios::binary);
  std::ofstream output(dest, std::ios::binary);
  if (!input || !output) throw Error("Cannot copy " + source + " to " + dest);
  char buffer[65536];
  while (input.read(buffer, sizeof(buffer)) || input.gcount()) output.write(buffer, input.gcount());
  output.close();
  if (input.bad() || !output) throw Error("Copy failed: " + source);
}
void atomic_text(const std::string &path, const std::string &text) {
  std::string temporary = path + ".tmp";
  try {
    std::ofstream output(temporary, std::ios::binary);
    output << text; output.close();
    if (!output) throw Error("Cannot write " + temporary);
    replace(temporary, path);
  } catch (...) { remove_file(temporary); throw; }
}
Json read_json(const std::string &path) {
  std::ifstream input(path);
  if (!input) throw Error("Cannot read " + path);
  Json result; input >> result; return result;
}
void atomic_json(const std::string &path, const Json &value) { atomic_text(path, value.dump(2) + "\n"); }
Json signature(const std::string &path) {
  struct stat info;
  if (stat(path.c_str(), &info)) throw Error(failure("Missing input", path));
#if defined(__APPLE__)
  int64_t ns = int64_t(info.st_mtimespec.tv_sec) * 1000000000 + info.st_mtimespec.tv_nsec;
#else
  int64_t ns = int64_t(info.st_mtim.tv_sec) * 1000000000 + info.st_mtim.tv_nsec;
#endif
  return {{"path", absolute(path)}, {"size", info.st_size}, {"mtime_ns", ns}};
}
std::string executable_path(const char *argv0) {
  char buffer[PATH_MAX]; ssize_t length = readlink("/proc/self/exe", buffer, sizeof(buffer) - 1);
  if (length > 0) { buffer[length] = 0; return buffer; }
  char *resolved = realpath(argv0, nullptr);
  if (!resolved) throw Error("Cannot locate executable");
  std::string result = resolved; free(resolved); return result;
}
std::string temporary_directory(const std::string &prefix) {
  std::string pattern = prefix + "XXXXXX";
  std::vector<char> buffer(pattern.begin(), pattern.end()); buffer.push_back(0);
  char *result = mkdtemp(buffer.data());
  if (!result) throw Error(failure("Cannot create temporary directory", prefix));
  return result;
}
void install_signals() {
  signal(SIGINT, signal_handler); signal(SIGTERM, signal_handler); signal(SIGPIPE, SIG_IGN);
}
void run(const std::string &executable, const std::vector<std::string> &arguments,
         const std::function<void(const std::string &)> &stderr_line, const std::string &stdout_path) {
  if (interrupted) throw Error("Interrupted");
  int pipes[2];
  if (pipe(pipes)) throw Error("Cannot create subprocess pipe");
  std::cout.flush(); std::cerr.flush(); fflush(nullptr);
  pid_t pid = fork();
  if (pid < 0) { close(pipes[0]); close(pipes[1]); throw Error("Cannot fork subprocess"); }
  if (pid == 0) {
    setpgid(0, 0);
    signal(SIGINT, SIG_DFL); signal(SIGTERM, SIG_DFL); signal(SIGPIPE, SIG_DFL);
    dup2(pipes[1], STDERR_FILENO); close(pipes[0]); close(pipes[1]);
    if (!stdout_path.empty()) {
      int output = open(stdout_path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0666);
      if (output < 0) _exit(126);
      dup2(output, STDOUT_FILENO); close(output);
    }
    setenv("MEGAHIT_HUMAN_READABLE_OUTPUT", "0", 1);
    std::vector<char *> argv; argv.push_back(const_cast<char *>(executable.c_str()));
    for (const auto &argument : arguments) argv.push_back(const_cast<char *>(argument.c_str()));
    argv.push_back(nullptr);
    execv(executable.c_str(), argv.data());
    dprintf(STDERR_FILENO, "Cannot execute %s: %s\n", executable.c_str(), std::strerror(errno));
    _exit(127);
  }
  setpgid(pid, pid); active_child = pid; close(pipes[1]);
  if (interrupted) kill(-pid, SIGTERM);
  FILE *stream = fdopen(pipes[0], "r");
  try {
    if (!stream) throw Error("Cannot read subprocess output");
    char *line = nullptr; size_t capacity = 0; ssize_t length;
    while ((length = getline(&line, &capacity, stream)) >= 0) {
      while (length && (line[length - 1] == '\n' || line[length - 1] == '\r')) --length;
      stderr_line(std::string(line, length));
    }
    free(line); fclose(stream); stream = nullptr;
    int result = wait_child(pid); active_child = -1;
    if (result != 0 || interrupted) throw Error("Subcommand failed (exit " + std::to_string(result) + "): " +
                                               (arguments.empty() ? executable : arguments[0]));
  } catch (...) {
    if (active_child == pid) { kill(-pid, SIGTERM); wait_child(pid); active_child = -1; }
    if (stream) fclose(stream);
    throw;
  }
}
void decompress(const std::string &source, const std::string &output) {
  int destination = open(output.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0666);
  if (destination < 0) throw Error(failure("Cannot open decompression output", output));
  try {
    char buffer[65536];
    if (ends_with(source, ".bz2")) {
      FILE *file = fopen(source.c_str(), "rb");
      if (!file) throw Error("Cannot open bzip2 input " + source);
      std::unique_ptr<FILE, int (*)(FILE *)> guard(file, fclose);
      std::vector<char> unused;
      while (true) {
        int error;
        BZFILE *bz = BZ2_bzReadOpen(&error, file, 0, 0,
                                   unused.empty() ? nullptr : unused.data(), unused.size());
        if (error != BZ_OK) throw Error("Cannot open bzip2 stream: " + source);
        try {
          int count;
          do {
            count = BZ2_bzRead(&error, bz, buffer, sizeof(buffer));
            if (error != BZ_OK && error != BZ_STREAM_END) throw Error("Corrupt bzip2 input: " + source);
            if (count) write_all(destination, buffer, count);
          } while (error == BZ_OK);
          void *remaining; int count_unused;
          BZ2_bzReadGetUnused(&error, bz, &remaining, &count_unused);
          if (error != BZ_OK) throw Error("Cannot finish bzip2 stream: " + source);
          std::vector<char> next;
          if (count_unused) next.assign(static_cast<char *>(remaining), static_cast<char *>(remaining) + count_unused);
          unused.swap(next);
        } catch (...) { BZ2_bzReadClose(&error, bz); throw; }
        BZ2_bzReadClose(&error, bz);
        if (unused.empty()) {
          int next = fgetc(file);
          if (next == EOF) break;
          ungetc(next, file);
        }
      }
    } else {
      gzFile input = gzopen(source.c_str(), "rb");
      if (!input) throw Error("Cannot open compressed input " + source);
      try {
        int count;
        while ((count = gzread(input, buffer, sizeof(buffer))) > 0) write_all(destination, buffer, count);
        if (count < 0) throw Error("Corrupt gzip input: " + source);
      } catch (...) { gzclose(input); throw; }
      if (gzclose(input) != Z_OK) throw Error("Cannot finish compressed input: " + source);
    }
    if (close(destination)) { destination = -1; throw Error("Cannot finish decompression output"); }
    destination = -1;
  } catch (...) { if (destination >= 0) close(destination); throw; }
}
void convert_fastq(const std::string &source, const std::string &output,
                   const std::string &pigz, int threads,
                   const std::function<void(const std::string &)> &log) {
  pid_t pid = -1; int input_fd = -1; gzFile gzip = nullptr;
  try {
    if (!pigz.empty() && threads > 1) {
      log("Merged FASTA compressor: bundled C pigz, " + std::to_string(threads) + " threads");
      int pipes[2]; if (pipe(pipes)) throw Error("Cannot create pigz pipe");
      pid = fork();
      if (pid < 0) { close(pipes[0]); close(pipes[1]); throw Error("Cannot fork pigz"); }
      if (pid == 0) {
        setpgid(0, 0); signal(SIGINT, SIG_DFL); signal(SIGTERM, SIG_DFL); signal(SIGPIPE, SIG_DFL);
        int destination = open(output.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0666);
        if (destination < 0) _exit(126);
        dup2(pipes[0], STDIN_FILENO); dup2(destination, STDOUT_FILENO);
        close(destination); close(pipes[0]); close(pipes[1]);
        std::string count = std::to_string(threads);
        execl(pigz.c_str(), pigz.c_str(), "-1", "-c", "-n", "-p", count.c_str(), nullptr);
        _exit(127);
      }
      setpgid(pid, pid); active_child = pid;
      close(pipes[0]); input_fd = pipes[1];
      fastq_records(source, [&](const std::string &data) { write_all(input_fd, data.data(), data.size()); });
      close(input_fd); input_fd = -1;
      int status = wait_child(pid); active_child = -1; pid = -1;
      if (status) throw Error("Bundled pigz failed (exit " + std::to_string(status) + ")");
    } else {
      log("Merged FASTA compressor: native C zlib (level 1)");
      gzip = gzopen(output.c_str(), "wb1");
      if (!gzip) throw Error("Cannot open gzip output " + output);
      gzbuffer(gzip, 1024 * 1024);
      fastq_records(source, [&](const std::string &data) {
        if (gzwrite(gzip, data.data(), data.size()) != static_cast<int>(data.size()))
          throw Error("Cannot write gzip output " + output);
      });
      int status = gzclose(gzip); gzip = nullptr;
      if (status != Z_OK) throw Error("Cannot finish gzip output " + output);
    }
    if (interrupted) throw Error("Interrupted");
  } catch (...) {
    if (input_fd >= 0) close(input_fd);
    if (gzip) gzclose(gzip);
    if (pid > 0) { kill(-pid, SIGTERM); wait_child(pid); active_child = -1; }
    remove_file(output); throw;
  }
}
void read_fasta(const std::string &path,
                const std::function<void(const std::string &, const std::string &)> &record) {
  std::ifstream input(path);
  if (!input) throw Error("Cannot open FASTA: " + path);
  std::string line, header, sequence;
  while (std::getline(input, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (!line.empty() && line[0] == '>') {
      if (!header.empty()) record(header, sequence);
      header = line; sequence.clear();
      if (line.find_first_not_of(" \t", 1) == std::string::npos) throw Error("Empty FASTA name: " + path);
    } else {
      auto first = line.find_first_not_of(" \t\r\n");
      if (first == std::string::npos) continue;
      if (header.empty()) throw Error("FASTA sequence without header: " + path);
      sequence += line.substr(first, line.find_last_not_of(" \t\r\n") - first + 1);
    }
  }
  if (input.bad()) throw Error("Cannot read FASTA: " + path);
  if (!header.empty()) record(header, sequence);
}
}  // namespace native
