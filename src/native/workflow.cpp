#include "platform.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sched.h>
#include <signal.h>
#include <sstream>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>
#include "definitions.h"
#include "minimap2/minimap.h"

namespace native {
namespace {
using Args = std::vector<std::string>;
template <typename T> std::string number(T value) { std::ostringstream out; out << value; return out.str(); }
std::vector<std::string> split(const std::string &value, char separator = ',') {
  std::vector<std::string> result; std::stringstream stream(value); std::string item;
  while (std::getline(stream, item, separator)) {
    if (item.empty()) throw Error("Empty item in list: " + value);
    result.push_back(item);
  }
  if (result.empty() || value.back() == separator) throw Error("Empty list item: " + value);
  return result;
}
int integer(const std::string &value) {
  size_t used; long long result;
  try { result = std::stoll(value, &used); } catch (...) { throw Error("Invalid integer: " + value); }
  if (used != value.size() || result < std::numeric_limits<int>::min() || result > std::numeric_limits<int>::max())
    throw Error("Invalid integer: " + value);
  return static_cast<int>(result);
}
double real(const std::string &value) {
  size_t used; double result;
  try { result = std::stod(value, &used); } catch (...) { throw Error("Invalid number: " + value); }
  if (used != value.size() || !std::isfinite(result)) throw Error("Invalid number: " + value);
  return result;
}
int cpu_count() {
  cpu_set_t set;
  if (sched_getaffinity(0, sizeof(set), &set) == 0) return std::max(1, CPU_COUNT(&set));
  return std::max<long>(1, sysconf(_SC_NPROCESSORS_ONLN));
}
struct Options {
  std::string out_dir, temp_dir, out_prefix, presets;
  std::string merge_index = "hash";
  std::string merge_dedup = "sort";
  double merge_sort_memory = 128;
  double memory = 0.9, merge_similar = 0.95, prune_depth = 2,
         disconnect_ratio = 0.1, low_local_ratio = 0.2, flash_max_mismatch_density = 0.25;
  int min_contig_len = 200, k_min = 21, k_max = 141, k_step = 10, min_count = 2,
      max_tip_len = -1, bubble_level = 2, merge_len = 20, prune_level = 2,
      num_cpu_threads = 0, cleaning_rounds = 5, mem_flag = 1,
      flash_min_overlap = 10, flash_max_overlap = 65, flash_phred_offset = 33;
  bool test_mode = false, continue_mode = false, force_overwrite = false,
       auto_k = true, set_list_by_min_max_step = false, no_mercy = false,
       no_local = false, keep_tmp_files = false, kmin_1pass = false,
       debug_assembly = false, merge_pairs = true, verbose = false, hw_accel = true;
  std::vector<int> k_list{21, 29, 39, 59, 79, 99, 119, 141};
  std::vector<int> reduction_k_list, flanking_k_list;
  bool legacy_k_list = false;
  bool reduction_k_list_set = false, flanking_k_list_set = false;
  bool flanking_consensus_only = false;
  Args pe1, pe2, pe12, se;
  NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT(Options, out_dir, temp_dir, out_prefix, presets, merge_index,
      merge_dedup, merge_sort_memory,
      memory, merge_similar, prune_depth, disconnect_ratio, low_local_ratio,
      flash_max_mismatch_density, min_contig_len, k_min, k_max, k_step, min_count,
      max_tip_len, bubble_level, merge_len, prune_level, num_cpu_threads, cleaning_rounds,
      mem_flag, flash_min_overlap, flash_max_overlap, flash_phred_offset, test_mode,
      auto_k, set_list_by_min_max_step, no_mercy, no_local, keep_tmp_files, kmin_1pass,
      debug_assembly, merge_pairs, verbose, hw_accel, k_list, reduction_k_list,
      flanking_k_list, legacy_k_list, reduction_k_list_set, flanking_k_list_set,
      flanking_consensus_only, pe1, pe2, pe12, se)
};
const char *help = R"(MEGAHIT native C++ workflow (no Python runtime required)
Usage: megahit [options] {-1 mates1 -2 mates2 | --12 interleaved | -r reads}

Inputs (repeatable, comma-separated; FASTA/FASTQ, plain/gzip/bzip2):
  -1 FILES  -2 FILES            Paired-end libraries
  --12 FILES                   Interleaved paired-end libraries
  -r, --read FILES             Single-end libraries
  --merge-pairs                Merge paired FASTQ before binary-library creation [on]
  --no-merge-pairs             Use original mates without FLASH (also for paired FASTA)
  --flash-min-overlap N        Minimum overlap [10]
  --flash-max-overlap N        Overlap scoring length [65]
  --flash-max-mismatch-density X  Mismatch density [0.25]
  --flash-phred-offset N       Quality offset, 33 or 64 [33]

Assembly:
  --reduction-k-list LIST      [21,29,39,59]
  --flanking-k-list LIST       [21,31,41,61,81,101,141]
  --flanking-consensus-only    Flanking uses retained contigs + consensus, no read edges
  --k-list LIST                Legacy shared list: reduction k<=60, flanking all
  --k-min N --k-max N --k-step N  Override k-list; odd k, even step <=28
  --min-count N                Minimum edge multiplicity [2]
  --merge-index fm|hash        Exact contig-end overlap index [hash]
  --merge-dedup bloom|sort     Merge candidate redundancy filter [sort]
  --merge-sort-memory X       Sort record buffers in MiB [128]
  --min-contig-len N           Final minimum contig length [200]
  --no-mercy --no-local --kmin-1pass --no-hw-accel
  --bubble-level N             Bubble cleaning intensity [2]
  --merge-level LENGTH,SIM     Complex bubble merging [20,0.95]
  --prune-level N --prune-depth X  [2,2]
  --disconnect-ratio X --low-local-ratio X  [0.1,0.2]
  --max-tip-len N --cleaning-rounds N --debug-assembly
  --presets meta-sensitive|meta-large

Resources and output:
  -t, --num-cpu-threads N      [available logical CPUs]
  -m, --memory X               Bytes, or physical-memory fraction [0.9]
  --mem-flag N                 Graph builder memory mode [1]
  -o, --out-dir DIR            [./megahit_out]
  --out-prefix NAME --tmp-dir DIR --keep-tmp-files
  --continue                   Resume this native driver's checkpoints
  --force --verbose --test --no-human-readable-output
  -h, --help  -v, --version

Reduction/merge follows its own list, then consensus flanking restarts at the
first flanking k and follows its own list. Phase-specific options override the
shared list for that phase. Reduction steps must be <=28, flanking steps <=40.
The final executable contains
the assembly core, FLASH, Minimap2, and optional pigz; gzip fallback uses zlib.
Native checkpoints are separate from legacy Python checkpoints. Use a fresh
output directory when switching from the legacy driver.
)";

Options parse(int argc, char **argv) {
  Options opt;
  if (argc == 1) throw Error(help);
  for (int i = 1; i < argc; ++i) {
    std::string key = argv[i], attached;
    if (key.compare(0, 2, "--") == 0) {
      auto position = key.find('=');
      if (position != std::string::npos) { attached = key.substr(position + 1); key.resize(position); }
    } else if (key.size() > 2 && key[0] == '-') { attached = key.substr(2); key.resize(2); }
    auto value = [&]() {
      if (!attached.empty()) { std::string result = attached; attached.clear(); return result; }
      if (++i >= argc) throw Error("Missing value for " + key);
      return std::string(argv[i]);
    };
    auto inputs = [&](Args &dest) { for (const auto &path : split(value())) dest.push_back(absolute(path)); };
    if (key == "-h" || key == "--help") { std::cout << help; std::exit(0); }
    else if (key == "-v" || key == "--version") { std::cout << "MEGAHIT " PACKAGE_VERSION " (native)\n"; std::exit(0); }
    else if (key == "-1") inputs(opt.pe1);
    else if (key == "-2") inputs(opt.pe2);
    else if (key == "--12") inputs(opt.pe12);
    else if (key == "-r" || key == "--read") inputs(opt.se);
    else if (key == "-o" || key == "--out-dir") opt.out_dir = absolute(value());
    else if (key == "--tmp-dir") opt.temp_dir = absolute(value());
    else if (key == "--out-prefix") opt.out_prefix = value();
    else if (key == "--presets") opt.presets = value();
    else if (key == "--merge-index") {
      opt.merge_index = value();
      if (opt.merge_index != "fm" && opt.merge_index != "hash") throw Error("--merge-index must be fm or hash");
    }
    else if (key == "--merge-dedup") {
      opt.merge_dedup = value();
      if (opt.merge_dedup != "bloom" && opt.merge_dedup != "sort") throw Error("--merge-dedup must be bloom or sort");
    }
    else if (key == "--merge-sort-memory") {
      opt.merge_sort_memory = real(value());
      if (opt.merge_sort_memory < 0.0625 || opt.merge_sort_memory > 1048576)
        throw Error("--merge-sort-memory must be between 0.0625 and 1048576 MiB");
    }
    else if (key == "-t" || key == "--num-cpu-threads") opt.num_cpu_threads = integer(value());
    else if (key == "-m" || key == "--memory") opt.memory = real(value());
    else if (key == "--min-contig-len") opt.min_contig_len = integer(value());
    else if (key == "--min-count") opt.min_count = integer(value());
    else if (key == "--mem-flag") opt.mem_flag = integer(value());
    else if (key == "--flash-min-overlap") opt.flash_min_overlap = integer(value());
    else if (key == "--flash-max-overlap") opt.flash_max_overlap = integer(value());
    else if (key == "--flash-max-mismatch-density") opt.flash_max_mismatch_density = real(value());
    else if (key == "--flash-phred-offset") opt.flash_phred_offset = integer(value());
    else if (key == "--reduction-k-list" || key == "--flanking-k-list") {
      auto &list = key == "--reduction-k-list" ? opt.reduction_k_list : opt.flanking_k_list;
      if (key == "--reduction-k-list") opt.reduction_k_list_set = true;
      else opt.flanking_k_list_set = true;
      list.clear(); for (const auto &k : split(value())) list.push_back(integer(k));
      if (list.empty()) throw Error(key + " cannot be empty");
      opt.auto_k = false;
    }
    else if (key == "--k-list") {
      opt.legacy_k_list = true;
      opt.k_list.clear(); for (const auto &k : split(value())) opt.k_list.push_back(integer(k));
      std::sort(opt.k_list.begin(), opt.k_list.end()); opt.auto_k = false; opt.set_list_by_min_max_step = false;
    } else if (key == "--k-min" || key == "--k-max" || key == "--k-step") {
      opt.legacy_k_list = true;
      int k = integer(value());
      if (key == "--k-min") opt.k_min = k;
      else if (key == "--k-max") opt.k_max = k;
      else opt.k_step = k;
      opt.auto_k = false; opt.set_list_by_min_max_step = true;
    } else if (key == "--max-tip-len") opt.max_tip_len = integer(value());
    else if (key == "--bubble-level") opt.bubble_level = integer(value());
    else if (key == "--prune-level") opt.prune_level = integer(value());
    else if (key == "--prune-depth") opt.prune_depth = real(value());
    else if (key == "--disconnect-ratio") opt.disconnect_ratio = real(value());
    else if (key == "--low-local-ratio") opt.low_local_ratio = real(value());
    else if (key == "--cleaning-rounds") opt.cleaning_rounds = integer(value());
    else if (key == "--merge-level") {
      auto fields = split(value()); if (fields.size() != 2) throw Error("--merge-level requires LENGTH,SIMILARITY");
      opt.merge_len = static_cast<int>(real(fields[0])); opt.merge_similar = real(fields[1]);
    } else if (key == "--merge-pairs") opt.merge_pairs = true;
    else if (key == "--no-merge-pairs") opt.merge_pairs = false;
    else if (key == "--no-mercy") opt.no_mercy = true;
    else if (key == "--no-local") opt.no_local = true;
    else if (key == "--flanking-consensus-only") opt.flanking_consensus_only = true;
    else if (key == "--kmin-1pass") opt.kmin_1pass = true;
    else if (key == "--no-hw-accel") opt.hw_accel = false;
    else if (key == "--keep-tmp-files") opt.keep_tmp_files = true;
    else if (key == "--continue") opt.continue_mode = true;
    else if (key == "--force" || key == "-f") opt.force_overwrite = true;
    else if (key == "--verbose") opt.verbose = true;
    else if (key == "--test") opt.test_mode = true;
    else if (key == "--debug-assembly") opt.debug_assembly = true;
    else if (key == "--no-human-readable-output") {}
    else if (key == "--cpu-only" || key == "--no-low-local" || key == "--use-gpu")
      std::cerr << "Ignoring deprecated option " << key << '\n';
    else if (key == "--max-read-len" || key == "--gpu-mem" || key == "-l") {
      value(); std::cerr << "Ignoring deprecated option " << key << '\n';
    } else throw Error("Unknown option: " + key);
    if (!attached.empty()) throw Error("Unexpected value for " + key);
  }
  return opt;
}
struct EarlyTerminate { int k; };

class Driver {
  Options opt;
  std::string executable, runtime_directory, input_workspace, pigz;
  std::ofstream logfile;
  Json completed = Json::object();
  Json visited = Json::object();
  bool flanking = false, resume_valid = true;
  uint64_t host_mem = 0;
  std::vector<pid_t> producers;
  Args fifos;
  std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

  std::string phase() const { return flanking ? "flanking/" : ""; }
  std::string contig_dir() const { return opt.out_dir + "/intermediate_contigs/" + (flanking ? "flanking" : ""); }
  std::string contig(int k) const { return contig_dir() + (flanking ? "/" : "") + "k" + number(k); }
  std::string graph(int k) const {
    std::string directory = opt.temp_dir + "/" + phase() + "k" + number(k);
    mkdirs(directory); return directory + "/" + number(k);
  }
  std::string reduction_path(int k, const std::string &type) const {
    std::string directory = opt.temp_dir + "/k" + number(k) + "/" + type;
    mkdirs(directory); return directory + "/" + number(k);
  }
  std::string read_count(int k) const { return reduction_path(k, "read_count"); }
  int initial_k() const {
    return opt.reduction_k_list.empty() ? opt.flanking_k_list.front() : opt.reduction_k_list.front();
  }
  std::string current_edges(int k) const {
    for (const auto &type : {"merged_edges", "reduced_edges"}) {
      auto prefix = reduction_path(k, type);
      if (exists(prefix + ".edges.info")) return prefix;
    }
    return k == initial_k() ? graph(k) : read_count(k);
  }
  std::string lib() const { return opt.temp_dir + "/reads.lib"; }
  std::string ends(int k) const {
    return opt.out_dir + "/tmp/" + phase() + "k" + number(k) + "/k" + number(k) + ".contig_ends.fa";
  }
  std::string checkpoints() const { return opt.out_dir + "/checkpoints.native.json"; }
  std::string final_path() const {
    return opt.out_dir + "/" + (opt.out_prefix.empty() ? "final" : opt.out_prefix) + ".contigs.fa";
  }
  void log(const std::string &message, bool console = true) {
    std::time_t now = std::time(nullptr); std::tm stamp; localtime_r(&now, &stamp);
    char formatted[32]; std::strftime(formatted, sizeof(formatted), "%Y-%m-%d %H:%M:%S", &stamp);
    std::string line = std::string(formatted) + " - " + message;
    if (logfile) { logfile << line << '\n'; logfile.flush(); }
    if (console || opt.verbose) std::cerr << line << '\n';
  }
  void core(const Args &arguments, const std::string &message, bool verbose = false,
             const std::string &stdout_path = "") {
    log(message);
    Args command{"--core"}; command.insert(command.end(), arguments.begin(), arguments.end());
    log("command " + executable + " " + Json(command).dump(), false);
    logfile.flush();
    run(executable, command, [&](const std::string &line) { log(line, verbose); }, stdout_path);
  }
  void checkpoint(const std::string &name, const Args &outputs, const std::function<void()> &action) {
    std::string key = phase() + name;
    bool valid = opt.continue_mode && resume_valid && completed.value(key, false);
    for (const auto &path : outputs) if (!exists(path)) valid = false;
    if (valid) { visited[key] = true; log("Reuse native checkpoint: " + key); return; }
    if (opt.continue_mode) resume_valid = false;
    // Invalidate the old tail before work starts, so a failed rerun cannot
    // mistake partial files for outputs of an older completed checkpoint.
    completed = visited;
    atomic_json(checkpoints(), completed);
    action();
    for (const auto &path : outputs) if (!exists(path)) throw Error("Stage did not produce " + path);
    visited[key] = true;
    completed = visited;
    atomic_json(checkpoints(), completed);
  }
  Args graph_options(int k, const std::string &output, int from = 0) const {
    return {"--host_mem", number(host_mem), "--mem_flag", number(opt.mem_flag),
            "--output_prefix", output, "--num_cpu_threads", number(opt.num_cpu_threads),
            "-k", number(k), "--kmer_from", number(from)};
  }
  void stop_producers() {
    for (pid_t pid : producers) {
      if (pid <= 0) continue;
      int status;
      if (waitpid(pid, &status, WNOHANG) == 0) {
        kill(pid, SIGTERM);
        while (waitpid(pid, &status, 0) < 0 && errno == EINTR) {}
      }
    }
    producers.clear();
    for (const auto &fifo : fifos) remove_file(fifo);
    fifos.clear();
    if (!input_workspace.empty()) { remove_tree(input_workspace); input_workspace.clear(); }
  }
  void validate() {
    if (!opt.presets.empty()) {
      opt.auto_k = !opt.reduction_k_list_set && !opt.flanking_k_list_set;
      opt.legacy_k_list = true; opt.min_count = 1;
      if (opt.presets == "meta-sensitive") {
        opt.k_list = {21,29,39,49,59,69,79,89,99,109,119,129,141}; opt.set_list_by_min_max_step = false;
      } else if (opt.presets == "meta-large") {
        opt.k_min = 27; opt.k_max = 127; opt.k_step = 10; opt.set_list_by_min_max_step = true;
      } else throw Error("Unknown preset: " + opt.presets);
    }
    if (opt.set_list_by_min_max_step) {
      if (opt.k_step <= 0 || opt.k_step > 28 || opt.k_step % 2 || opt.k_min > opt.k_max ||
          opt.k_min < 15 || opt.k_max > static_cast<int>(kMaxK) || opt.k_min % 2 == 0 || opt.k_max % 2 == 0)
        throw Error("k-step must be positive/even and <=28, with k-min <= k-max");
      opt.k_list.clear();
      for (int k = opt.k_min; k < opt.k_max; k += opt.k_step) opt.k_list.push_back(k);
      opt.k_list.push_back(opt.k_max);
    }
    if (opt.k_list.empty()) throw Error("k-list cannot be empty");
    for (size_t i = 0; i < opt.k_list.size(); ++i) {
      int k = opt.k_list[i];
      if (k < 15 || k > static_cast<int>(kMaxK) || k % 2 == 0) throw Error("All k values must be odd and in [15,255]");
      if (i && (k <= opt.k_list[i - 1] || k - opt.k_list[i - 1] > 28))
        throw Error("k-list must increase strictly, with adjacent difference <=28");
    }
    opt.k_min = opt.k_list.front(); opt.k_max = opt.k_list.back();
    if (opt.reduction_k_list.empty()) {
      if (opt.legacy_k_list) {
        for (int k : opt.k_list) if (k <= 60) opt.reduction_k_list.push_back(k);
      } else opt.reduction_k_list = {21,29,39,59};
    }
    if (opt.flanking_k_list.empty())
      opt.flanking_k_list = opt.legacy_k_list ? opt.k_list : std::vector<int>{21,31,41,61,81,101,141};
    auto check_list = [&](const std::vector<int> &list, const std::string &name, int max_step) {
      for (size_t i = 0; i < list.size(); ++i) {
        int k = list[i];
        if (k < 15 || k > static_cast<int>(kMaxK) || k % 2 == 0)
          throw Error(name + " values must be odd and in [15,255]");
        if (i && (k <= list[i - 1] || k - list[i - 1] > max_step))
          throw Error(name + " must increase strictly, with adjacent difference <=" + number(max_step));
      }
    };
    check_list(opt.reduction_k_list, "reduction-k-list", 28);
    check_list(opt.flanking_k_list, "flanking-k-list", 40);
    if (opt.min_count <= 0 || opt.min_contig_len < 0) throw Error("Invalid minimum count/contig length");
    if (opt.min_count == 1) { opt.no_mercy = true; opt.kmin_1pass = true; }
    if (opt.memory <= 0 || !std::isfinite(opt.memory)) throw Error("Memory must be positive");
    if (opt.memory < 1) {
      long pages = sysconf(_SC_PHYS_PAGES), page_size = sysconf(_SC_PAGESIZE);
      if (pages <= 0 || page_size <= 0) throw Error("Cannot detect physical memory; pass a byte count with -m");
      host_mem = static_cast<uint64_t>(std::floor(opt.memory * pages * page_size));
    } else {
      if (opt.memory > static_cast<double>(std::numeric_limits<int64_t>::max())) throw Error("Memory byte count too large");
      host_mem = static_cast<uint64_t>(opt.memory);
    }
    if (!host_mem) throw Error("Memory byte count must be at least one");
    if (opt.num_cpu_threads < 0) throw Error("Thread count cannot be negative");
    if (opt.num_cpu_threads == 0 || opt.num_cpu_threads > cpu_count()) opt.num_cpu_threads = cpu_count();
    if (opt.prune_level < 0 || opt.prune_level > 3 || opt.merge_len < 0 ||
        opt.merge_similar < 0 || opt.merge_similar > 1 || opt.disconnect_ratio < 0 ||
        opt.disconnect_ratio > .5 || opt.low_local_ratio <= 0 || opt.low_local_ratio > .5 ||
        opt.cleaning_rounds <= 0) throw Error("Invalid assembly cleaning parameters");
    opt.bubble_level = std::max(0, std::min(2, opt.bubble_level));
    if (opt.prune_depth < 0 && opt.prune_level < 3) opt.prune_depth = opt.min_count;
    if (opt.flash_min_overlap < 1 || opt.flash_max_overlap < opt.flash_min_overlap ||
        opt.flash_max_mismatch_density < 0 || opt.flash_max_mismatch_density > 1 ||
        (opt.flash_phred_offset != 33 && opt.flash_phred_offset != 64)) throw Error("Invalid FLASH parameters");
    if (opt.out_prefix.find('/') != std::string::npos) throw Error("out-prefix must be a file name, not a path");
  }
  Json inputs_signature() const {
    Json result = Json::array();
    for (const auto *files : {&opt.pe1, &opt.pe2, &opt.pe12, &opt.se})
      for (const auto &path : *files) {
        if (path.find_first_of("\r\n") != std::string::npos) throw Error("Input paths cannot contain line breaks");
        result.push_back(signature(path));
      }
    return result;
  }
  void setup() {
    if (opt.out_dir.empty()) opt.out_dir = absolute("megahit_out");
    if (opt.out_dir == "/") throw Error("The filesystem root cannot be an output directory");
    if (opt.continue_mode) {
      Json saved = read_json(opt.out_dir + "/options.json");
      if (saved.value("workflow_driver", std::string()) != "native-cpp-v1")
        throw Error("Legacy Python checkpoints cannot be resumed by the native driver; use a fresh output directory");
      std::string output = opt.out_dir;
      opt = saved.get<Options>(); opt.out_dir = output; opt.continue_mode = true;
      if (!saved.contains("reduction_k_list") && !saved.contains("flanking_k_list"))
        opt.legacy_k_list = true;
      // Native runs predating this option used Bloom deduplication.
      if (!saved.contains("merge_dedup")) opt.merge_dedup = "bloom";
      if (exists(checkpoints())) completed = read_json(checkpoints());
    } else {
      if (exists(opt.out_dir) && !opt.force_overwrite) throw Error("Output directory already exists: " + opt.out_dir);
      mkdirs(opt.out_dir);
      if (opt.temp_dir.empty()) opt.temp_dir = opt.out_dir + "/tmp";
      else { mkdirs(opt.temp_dir); opt.temp_dir = temporary_directory(opt.temp_dir + "/megahit_tmp_"); }
    }
    validate();
    if (opt.test_mode) {
      opt.merge_pairs = false;  // The bundled test libraries are FASTA.
      Args candidates{parent(executable) + "/../test_data", parent(executable) + "/../share/megahit/test_data",
                      MEGAHIT_TEST_DATA_DIR};
      std::string data;
      for (const auto &candidate : candidates) if (exists(candidate + "/r1.il.fa.gz")) { data = absolute(candidate); break; }
      if (data.empty()) throw Error("Test data not found next to the executable or installation");
      opt.pe12 = {data + "/r1.il.fa.gz"};
      if (exists(data + "/r2.il.fa.bz2")) opt.pe12.push_back(data + "/r2.il.fa.bz2");
      opt.pe1 = {data + "/r3_1.fa"}; opt.pe2 = {data + "/r3_2.fa"};
      opt.se = {data + "/r4.fa", data + "/loop.fa"};
    }
    if (opt.pe1.size() != opt.pe2.size()) throw Error("Number of paired-end files does not match");
    if (opt.pe1.empty() && opt.pe12.empty() && opt.se.empty()) throw Error("No read libraries provided");
    Json fingerprint = inputs_signature();
    if (opt.continue_mode && read_json(opt.out_dir + "/options.json").value("input_signatures", Json()) != fingerprint)
      throw Error("Input files changed since the native run; use a fresh output directory");
    if (!opt.continue_mode) {
      remove_file(checkpoints());
      remove_file(opt.out_dir + "/done");
    }
    mkdirs(opt.temp_dir); mkdirs(contig_dir()); mkdirs(opt.out_dir + "/tmp");
    logfile.open(opt.out_dir + "/" + (opt.out_prefix.empty() ? "log" : opt.out_prefix + ".log"), std::ios::app);
    if (!logfile) throw Error("Cannot open assembly log");
    Json saved = opt; saved["workflow_driver"] = "native-cpp-v1"; saved["input_signatures"] = fingerprint;
    atomic_json(opt.out_dir + "/options.json", saved);
    log("MEGAHIT " PACKAGE_VERSION " native C++ workflow; no Python runtime");
    log("Threads: " + number(opt.num_cpu_threads) + "; memory bytes: " + number(host_mem));
    if (!opt.hw_accel) log("Generic CPU path requested; native kernels use their existing runtime dispatch");
  }
  std::string input_path(const std::string &path, const std::string &type, size_t index) const {
    return ends_with(path, ".gz") || ends_with(path, ".bz2") || path.find_first_of(" \t") != std::string::npos ?
        input_workspace + "/inpipe." + type + "." + number(index) : path;
  }
  void library() {
    checkpoint("read-library", {lib() + ".bin", lib() + ".lib_info"}, [&]() {
      input_workspace = temporary_directory("/tmp/megahit-input-");
      std::ostringstream descriptor;
      auto add = [&](const std::string &path, const std::string &type, size_t index) {
        std::string converted = input_path(path, type, index);
        if (converted != path) {
          if (!ends_with(path, ".gz") && !ends_with(path, ".bz2")) {
            if (symlink(path.c_str(), converted.c_str())) throw Error("Cannot create read input alias");
            return converted;
          }
          remove_file(converted);
          if (mkfifo(converted.c_str(), 0600)) throw Error("Cannot create input FIFO " + converted);
          fifos.push_back(converted);
          logfile.flush();
          pid_t pid = fork();
          if (pid < 0) throw Error("Cannot fork input decompressor");
          if (pid == 0) {
            signal(SIGTERM, SIG_DFL); signal(SIGINT, SIG_DFL);
            execl(executable.c_str(), executable.c_str(), "--native-decompress", path.c_str(), converted.c_str(), nullptr);
            _exit(127);
          }
          producers.push_back(pid);
        }
        return converted;
      };
      try {
        if (opt.merge_pairs && (!opt.pe1.empty() || !opt.pe12.empty())) {
          Args inputs = read_json(opt.out_dir + "/tmp/flash/assembly_inputs.json").get<Args>();
          if (inputs.empty() && opt.se.empty()) throw Error("FLASH produced no reads for the assembly library");
          for (size_t i = 0; i < inputs.size(); ++i)
            descriptor << inputs[i] << "\nse " << add(inputs[i], "flash", i) << '\n';
        } else {
          for (size_t i = 0; i < opt.pe12.size(); ++i)
            descriptor << opt.pe12[i] << "\ninterleaved " << add(opt.pe12[i], "pe12", i) << '\n';
          for (size_t i = 0; i < opt.pe1.size(); ++i)
            descriptor << opt.pe1[i] << ',' << opt.pe2[i] << "\npe " << add(opt.pe1[i], "pe1", i) << ' '
                       << add(opt.pe2[i], "pe2", i) << '\n';
        }
        for (size_t i = 0; i < opt.se.size(); ++i)
          descriptor << opt.se[i] << "\nse " << add(opt.se[i], "se", i) << '\n';
        atomic_text(lib(), descriptor.str());
        core({"buildlib", lib(), lib(), number(opt.num_cpu_threads)}, "Convert reads to binary library", true);
        for (pid_t &pid : producers) {
          int status; pid_t result;
          do { result = waitpid(pid, &status, 0); } while (result < 0 && errno == EINTR);
          pid = -1;
          if (result < 0 || !WIFEXITED(status) || WEXITSTATUS(status)) throw Error("Input decompression failed");
        }
        producers.clear();
        for (const auto &fifo : fifos) remove_file(fifo);
        fifos.clear();
        remove_tree(input_workspace); input_workspace.clear();
      } catch (...) { stop_producers(); throw; }
    });
    if (opt.auto_k && opt.legacy_k_list && opt.k_list.size() > 1) {
      std::ifstream info(lib() + ".lib_info"); std::string line; int row = 0, maximum = 0;
      while (std::getline(info, line)) {
        if (row >= 2 && row % 2 == 0) { std::istringstream values(line); int64_t begin, end; int length;
          if (values >> begin >> end >> length) maximum = std::max(maximum, length); }
        ++row;
      }
      std::vector<int> adjusted;
      for (int k : opt.k_list) if (k < maximum + 20) adjusted.push_back(k);
      if (!adjusted.empty()) {
        opt.k_list = adjusted; opt.k_min = adjusted.front(); opt.k_max = adjusted.back();
        opt.flanking_k_list = adjusted; opt.reduction_k_list.clear();
        for (int k : adjusted) if (k <= 60) opt.reduction_k_list.push_back(k);
      }
    }
  }
  std::string compressor() {
    if (!pigz.empty() || !bundled_pigz_size) return pigz;
    try {
      const char *temporary = std::getenv("TMPDIR");
      runtime_directory = temporary_directory(std::string(temporary ? temporary : "/tmp") + "/megahit-native-");
      pigz = runtime_directory + "/pigz";
      std::ofstream output(pigz, std::ios::binary);
      output.write(reinterpret_cast<const char *>(bundled_pigz), bundled_pigz_size); output.close();
      if (!output || chmod(pigz.c_str(), 0700)) throw Error("Cannot unpack native pigz");
    } catch (const std::exception &error) { log(std::string("Optional pigz unavailable: ") + error.what()); pigz.clear(); }
    return pigz;
  }
  void merge_reads() {
    if (!opt.merge_pairs || (opt.pe1.empty() && opt.pe12.empty())) {
      if (!opt.continue_mode) {
        remove_file(opt.out_dir + "/tmp/flash/merged.json");
        remove_file(opt.out_dir + "/tmp/flash/assembly_inputs.json");
      }
      return;
    }
    std::string directory = opt.out_dir + "/tmp/flash"; mkdirs(directory);
    std::string manifest = directory + "/merged.json";
    std::string assembly_manifest = directory + "/assembly_inputs.json";
    size_t libraries = opt.pe1.size() + opt.pe12.size();
    Args outputs;
    for (size_t i = 0; i < libraries; ++i) outputs.push_back(directory + "/library." + number(i) + ".merged.fa.gz");
    Args expected = outputs; expected.push_back(manifest); expected.push_back(assembly_manifest);
    for (size_t i = 0; i < libraries; ++i)
      for (int mate = 1; mate <= 2; ++mate)
        expected.push_back(directory + "/library." + number(i) + ".unmerged." + number(mate) + ".fa.gz");
    checkpoint("flash-assembly-inputs-v2", expected, [&]() {
      remove_file(manifest);
      remove_file(assembly_manifest);
      Args assembly_inputs;
      for (size_t i = 0; i < libraries; ++i) {
        std::string work = temporary_directory(opt.temp_dir + "/flash-");
        try {
          bool interleaved = i >= opt.pe1.size();
          Args inputs = interleaved ? Args{opt.pe12[i - opt.pe1.size()]} : Args{opt.pe1[i], opt.pe2[i]};
          for (size_t j = 0; j < inputs.size(); ++j) if (ends_with(inputs[j], ".bz2")) {
            std::string unpacked = work + "/input." + number(j) + ".fastq";
            decompress(inputs[j], unpacked); inputs[j] = unpacked;
          }
          std::string prefix = "library." + number(i);
          Args args{"flash", "-t", number(opt.num_cpu_threads), "-m", number(opt.flash_min_overlap),
                    "-M", number(opt.flash_max_overlap), "-x", number(opt.flash_max_mismatch_density),
                    "-p", number(opt.flash_phred_offset), "-d", work, "-o", prefix};
          if (interleaved) args.push_back("--interleaved-input");
          args.insert(args.end(), inputs.begin(), inputs.end());
          core(args, "Merge paired-end library " + number(i) + " with FLASH", true);
          auto begin = std::chrono::steady_clock::now();
          auto export_reads = [&](const std::string &source, const std::string &target) {
            try {
              convert_fastq(source, target + ".tmp", opt.num_cpu_threads > 1 ? compressor() : "",
                            opt.num_cpu_threads, [&](const std::string &message) { log(message); });
              replace(target + ".tmp", target);
            } catch (...) { remove_file(target + ".tmp"); throw; }
            // Empty gzip files are valid exports, but not binary-library inputs.
            if (size(source)) assembly_inputs.push_back(target);
          };
          export_reads(work + "/" + prefix + ".extendedFrags.fastq", outputs[i]);
          log("Merged FASTA: " + outputs[i] + " (conversion " + number(std::chrono::duration<double>(
              std::chrono::steady_clock::now() - begin).count()) + " seconds)");
          for (int mate = 1; mate <= 2; ++mate)
            export_reads(work + "/" + prefix + ".notCombined_" + number(mate) + ".fastq",
                         directory + "/" + prefix + ".unmerged." + number(mate) + ".fa.gz");
          remove_tree(work);
        } catch (...) { remove_file(outputs[i] + ".tmp"); remove_tree(work); throw; }
      }
      atomic_json(manifest, outputs);
      atomic_json(assembly_manifest, assembly_inputs);
    });
  }
  void count(int k, bool first) {
    std::string output = first ? graph(k) : read_count(k);
    checkpoint("count-k" + number(k), {output + ".edges.info"}, [&]() {
      core({"count", "-k", number(k), "-m", number(opt.min_count), "--host_mem", number(host_mem),
            "--mem_flag", number(opt.mem_flag), "--output_prefix", output,
            "--num_cpu_threads", number(opt.num_cpu_threads), "--read_lib_file", lib()},
           "Extract read-counted (k+1)-mers for k = " + number(k));
    });
  }
  void graph_edges(int k, const std::string &input, const std::string &output, const std::string &name) {
    checkpoint(name, {output + ".sdbg_info"}, [&]() {
      Args args{"seq2sdbg"}; auto common = graph_options(k, output);
      args.insert(args.end(), common.begin(), common.end());
      args.insert(args.end(), {"--input_prefix", input});
      core(args, "Build graph for k = " + number(k) + " from " + input);
    });
  }
  void clear_edges(const std::string &prefix) {
    if (!directory(parent(prefix))) return;
    std::string start = base(prefix) + ".edges.";
    for (const auto &name : entries(parent(prefix))) if (name.compare(0, start.size(), start) == 0)
      remove_file(parent(prefix) + "/" + name);
  }
  std::string reduce(int from, int to) {
    std::string input = from == initial_k() ? graph(from) : reduction_path(from, "merged_edges");
    if (!exists(input + ".edges.info") && from != initial_k()) input = read_count(from);
    std::string longer = read_count(to), output = reduction_path(from, "reduced_edges");
    checkpoint("reduce-k" + number(from) + "-k" + number(to), {output + ".edges.info"}, [&]() {
      if (!exists(input + ".edges.info") || !exists(longer + ".edges.info")) throw Error("Missing reduction input edges");
      clear_edges(output);
      core({"prefixbf", "--mode", "reduce", "--short_prefix", input, "--long_prefix", longer,
            "--output_prefix", output, "--short_k", number(from), "--long_k", number(to),
            "--num_cpu_threads", number(opt.num_cpu_threads), "--human_readable", "0",
            "--max_bloom_gb", "60", "--reserve_gb", "1", "--host_mem", number(host_mem)},
           "Reduce edges for k = " + number(from) + " using k = " + number(to), true);
    });
    return output;
  }
  void split_contigs(int from, int to) {
    std::string prefix = contig(from);
    std::string shorter = prefix + ".lt_k" + number(to) + "_plus1.contigs.fa";
    std::string longer = prefix + ".ge_k" + number(to) + "_plus1.contigs.fa";
    std::ofstream short_out(shorter), long_out(longer);
    uint64_t short_count = 0, short_bases = 0, long_count = 0, long_bases = 0;
    read_fasta(prefix + ".contigs.fa", [&](const std::string &header, const std::string &sequence) {
      bool small = sequence.size() < static_cast<size_t>(to + 1);
      std::ostream &dest = small ? short_out : long_out;
      dest << header << '\n' << sequence << '\n';
      if (small) { ++short_count; short_bases += sequence.size(); }
      else { ++long_count; long_bases += sequence.size(); }
    });
    short_out.close(); long_out.close();
    if (!short_out || !long_out) throw Error("Cannot split intermediate contigs");
    log("Split k = " + number(from) + " contigs at length " + number(to + 1) + ": short " +
        number(short_count) + " contigs/" + number(short_bases) + " bp, long " + number(long_count) +
        " contigs/" + number(long_bases) + " bp");
  }
  void merge_edges(int from, int to) {
    std::string prefix = contig(from), output = reduction_path(to, "merged_edges");
    checkpoint("merge-k" + number(from) + "-k" + number(to), {output + ".edges.info"}, [&]() {
      clear_edges(output);
      core({"prefixbf", "--mode", "merge", "--long_prefix", read_count(to), "--output_prefix", output,
            "--contig", prefix + ".lt_k" + number(to) + "_plus1.contigs.fa",
            "--long_contig", prefix + ".ge_k" + number(to) + "_plus1.contigs.fa",
            "--short_k", number(from), "--long_k", number(to), "--num_cpu_threads",
            number(opt.num_cpu_threads), "--human_readable", "0", "--merge_index", opt.merge_index,
            "--merge_dedup", opt.merge_dedup, "--sort_memory_mb", number(opt.merge_sort_memory),
            "--host_mem", number(host_mem)},
           "Merge k = " + number(from) + " assembly edges into k = " + number(to), true);
    });
  }
  void assemble(int k, const std::string &input, int next = 0) {
    int final_k = opt.flanking_k_list.back();
    if (!flanking && !opt.reduction_k_list.empty()) final_k = std::max(final_k, opt.reduction_k_list.back());
    Args outputs{contig(k) + ".contigs.fa", contig(k) + ".contigs.fa.info",
                 contig(k) + ".bubble_seq.fa", contig(k) + ".bubble_seq.fa.info"};
    if (next) { outputs.push_back(contig(k) + ".lt_k" + number(next) + "_plus1.contigs.fa");
                outputs.push_back(contig(k) + ".ge_k" + number(next) + "_plus1.contigs.fa"); }
    checkpoint("assemble-k" + number(k), outputs, [&]() {
      mkdirs(contig_dir());
      remove_file(contig(k) + ".addi.fa");
      remove_file(contig(k) + ".addi.fa.info");
      int standalone = std::max(std::min(final_k * 3 - 1, int(opt.min_contig_len * 1.5)), opt.min_contig_len);
      if (opt.max_tip_len >= 0) standalone = std::max(opt.max_tip_len + final_k - 1, opt.min_contig_len);
      double tip = opt.max_tip_len;
      if (tip == -1 && k * 3 - 1 > opt.min_contig_len * 1.5) tip = std::max(1.0, opt.min_contig_len * 1.5 + 1 - k);
      Args args{"assemble", "-s", input, "-o", contig(k), "-t", number(opt.num_cpu_threads),
                "--min_standalone", number(next ? 0 : standalone), "--prune_level", number(next ? 0 : opt.prune_level),
                "--merge_len", number(next ? 0 : opt.merge_len), "--merge_similar", number(opt.merge_similar),
                "--cleaning_rounds", number(next ? 1 : opt.cleaning_rounds),
                "--disconnect_ratio", number(opt.disconnect_ratio), "--low_local_ratio", number(opt.low_local_ratio),
                "--min_depth", number(next ? 0 : opt.prune_depth), "--bubble_level", number(next ? 0 : opt.bubble_level),
                "--max_tip_len", number(next ? 0 : tip)};
      if (k < final_k && !next) args.push_back("--careful_bubble");
      if (k == final_k) args.push_back("--is_final_round");
      if (opt.no_local && (flanking || opt.reduction_k_list.empty())) args.push_back("--output_standalone");
      core(args, "Assemble contigs from SdBG for k = " + number(k));
      if (next) split_contigs(k, next);
    });
  }
  void graph_sequences(int k, int from) {
    std::string output = graph(k);
    checkpoint("sequence-graph-k" + number(k), {output + ".sdbg_info"}, [&]() {
      Args args{"seq2sdbg"}; auto common = graph_options(k, output, from);
      args.insert(args.end(), common.begin(), common.end());
      uint64_t bases = 0;
      bool consensus_only = opt.flanking_consensus_only && from;
      if (!consensus_only && exists(output + ".edges.0")) {
        args.insert(args.end(), {"--input_prefix", output});
        for (int shard = 0; exists(output + ".edges." + number(shard)); ++shard)
          bases += size(output + ".edges." + number(shard));
      }
      if (from) {
        std::string prefix = contig(from);
        for (const auto &input : std::vector<std::pair<std::string, std::string>>{
             {"--addi_contig", ".addi.fa"}, {"--local_contig", ".local.fa"}}) {
          if (input.first == "--local_contig" && opt.no_local) continue;
          if (exists(prefix + input.second)) { args.insert(args.end(), {input.first, prefix + input.second}); bases += size(prefix + input.second); }
        }
        if (exists(prefix + ".contigs.fa")) {
          args.insert(args.end(), {"--contig", prefix + ".contigs.fa", "--bubble", prefix + ".bubble_seq.fa"});
          if (consensus_only) {
            bases += size(prefix + ".contigs.fa");
            if (exists(prefix + ".bubble_seq.fa")) bases += size(prefix + ".bubble_seq.fa");
          }
        }
      }
      if (from && !bases) throw EarlyTerminate{from};
      if (!opt.no_mercy && k == opt.flanking_k_list.front()) args.push_back("--need_mercy");
      core(args, "Build graph for k = " + number(k));
    });
  }
  void extract_ends(int k) {
    std::string output = ends(k); mkdirs(parent(output));
    std::ostringstream fasta;
    uint64_t whole = 0, divided = 0;
    read_fasta(contig(k) + ".contigs.fa", [&](const std::string &header, const std::string &sequence) {
      if (sequence.size() < 580) { fasta << header << '\n' << sequence << '\n'; ++whole; }
      else {
        size_t begin = header.find_first_not_of(" \t", 1);
        std::string name = header.substr(begin, header.find_first_of(" \t", begin) - begin);
        fasta << '>' << name << "|start source_length=" << sequence.size() << " region=1-290\n" << sequence.substr(0, 290) << '\n';
        fasta << '>' << name << "|end source_length=" << sequence.size() << " region=" << sequence.size() - 289 << '-'
              << sequence.size() << '\n' << sequence.substr(sequence.size() - 290) << '\n'; ++divided;
      }
    });
    bool same = false;
    if (exists(output)) { std::ifstream previous(output); std::ostringstream data; data << previous.rdbuf(); same = data.str() == fasta.str(); }
    if (!same) atomic_text(output, fasta.str());
    log("Contig ends from k = " + number(k) + ": " + number(whole) + " whole, " + number(divided) + " split");
  }
  void local(int k, int to) {
    std::string output = contig(k) + ".local.fa";
    checkpoint("consensus-local-read-support-primary-580-all-flash-v8-k" + number(k) + "-k" + number(to), {output, output + ".info"}, [&]() {
      log("Consensus local contigs for k = " + number(k) + " -> " + number(to));
      std::string manifest = opt.out_dir + "/tmp/flash/assembly_inputs.json";
      if (!opt.merge_pairs || !exists(manifest)) {
        atomic_text(output, ""); atomic_text(output + ".info", "0 0\n");
        log("No FLASH assembly read libraries; write empty consensus local contigs"); return;
      }
      Args queries = read_json(manifest).get<Args>();
      extract_ends(k);
      if (queries.empty() || !size(ends(k))) { atomic_text(output, ""); atomic_text(output + ".info", "0 0\n"); return; }
      std::string mapping = parent(ends(k)) + "/flank_mapping"; mkdirs(mapping);
      Args params{"-x", "sr", "-a", "-Y", "--secondary=no",
                  "-t", number(opt.num_cpu_threads), "-I", "55G"};
      Json configuration{{"version", MM_VERSION}, {"parameters", params}, {"reference", signature(ends(k))}};
      configuration["queries"] = Json::array();
      Args alignments;
      for (size_t i = 0; i < queries.size(); ++i) {
        configuration["queries"].push_back(signature(queries[i]));
        alignments.push_back(mapping + "/library." + number(i) + ".sam");
      }
      bool cached = false;
      if (opt.continue_mode && exists(mapping + "/mapping.json")) {
        Json previous = read_json(mapping + "/mapping.json");
        cached = previous.value("configuration", Json()) == configuration && previous.value("outputs", Json()) == Json(alignments);
        for (const auto &alignment : alignments) if (!exists(alignment) || !size(alignment)) cached = false;
      }
      if (!cached) {
        remove_file(mapping + "/mapping.json");
        for (size_t i = 0; i < queries.size(); ++i) {
          Args args{"minimap2"}; args.insert(args.end(), params.begin(), params.end());
          args.insert(args.end(), {"-o", alignments[i] + ".tmp", ends(k), queries[i]});
          try { core(args, "Map assembly reads to k = " + number(k) + " contig ends", true); replace(alignments[i] + ".tmp", alignments[i]); }
          catch (...) { remove_file(alignments[i] + ".tmp"); throw; }
        }
        atomic_json(mapping + "/mapping.json", {{"configuration", configuration}, {"outputs", alignments}});
      }
      std::string extended = parent(ends(k)) + "/k" + number(k) + ".contig_ends.extended.fa";
      std::string stats = mapping + "/consensus.tsv";
      Args args{"endconsensus", ends(k), extended + ".tmp", stats + ".tmp"};
      args.insert(args.end(), alignments.begin(), alignments.end());
      try {
        core(args, "Build outward contig-end consensus at k = " + number(k), true);
        replace(extended + ".tmp", extended); replace(extended + ".tmp.info", extended + ".info"); replace(stats + ".tmp", stats);
      } catch (...) { remove_file(extended + ".tmp"); remove_file(extended + ".tmp.info"); remove_file(stats + ".tmp"); throw; }
      Json consensus{{"algorithm", "outward-plurality-read-support-primary-580-all-flash-v8"}, {"reference", signature(ends(k))}};
      consensus["alignments"] = Json::array();
      for (const auto &alignment : alignments) consensus["alignments"].push_back(signature(alignment));
      atomic_json(mapping + "/consensus.json", consensus);
      copy(extended, output + ".tmp"); copy(extended + ".info", output + ".info.tmp");
      replace(output + ".tmp", output); replace(output + ".info.tmp", output + ".info");
    });
  }
  void seed(int source_k, const std::string &source) {
    int k = opt.flanking_k_list.front(); std::string output = graph(k);
    checkpoint("seed-from-k" + number(source_k), {output + ".sdbg_info"}, [&]() {
      Args args{"seq2sdbg"}; auto common = graph_options(k, output, source_k);
      args.insert(args.end(), common.begin(), common.end());
      args.insert(args.end(), {"--contig", source + ".contigs.fa", "--bubble", source + ".bubble_seq.fa"});
      for (const auto &input : std::vector<std::pair<std::string, std::string>>{
           {"--addi_contig", ".addi.fa"}, {"--local_contig", ".local.fa"}}) {
        if (input.first == "--local_contig" && opt.no_local) continue;
        if (exists(source + input.second)) args.insert(args.end(), {input.first, source + input.second});
      }
      core(args, "Seed flanking k = " + number(k) + " from k = " + number(source_k) + " contigs");
    });
  }
  void iterate(int from, int to) {
    std::string output = graph(to);
    checkpoint("iterate-k" + number(from) + "-k" + number(to), {output + ".edges.info"}, [&]() {
      core({"iterate", "-c", contig(from) + ".contigs.fa", "-b", contig(from) + ".bubble_seq.fa",
            "-t", number(opt.num_cpu_threads), "-k", number(from), "-s", number(to - from),
            "-o", output, "-r", lib() + ".bin"},
           "Extract iterative edges from k = " + number(from) + " to " + number(to));
    });
  }
  void merge_final(int k) {
    checkpoint("final-contigs-k" + number(k), {final_path()}, [&]() {
      std::ofstream output(final_path() + ".tmp");
      std::vector<size_t> lengths;
      auto append = [&](const std::string &path) {
        read_fasta(path, [&](const std::string &header, const std::string &sequence) {
          if (sequence.size() >= static_cast<size_t>(opt.min_contig_len)) {
            output << header << '\n' << sequence << '\n'; lengths.push_back(sequence.size());
          }
        });
      };
      for (const auto &name : entries(contig_dir())) {
        for (int assembled_k : ((flanking || opt.reduction_k_list.empty()) ? opt.flanking_k_list : opt.reduction_k_list)) {
          if (assembled_k <= k && name == "k" + number(assembled_k) + ".final.contigs.fa")
            append(contig_dir() + "/" + name);
        }
      }
      append(contig(k) + ".contigs.fa");
      output.close(); if (!output) throw Error("Cannot write final contigs");
      replace(final_path() + ".tmp", final_path());
      uint64_t total = 0, cumulative = 0, n50 = 0;
      for (auto length : lengths) total += length;
      std::sort(lengths.rbegin(), lengths.rend());
      for (auto length : lengths) { cumulative += length; if (cumulative >= (total + 1) / 2) { n50 = length; break; } }
      log(number(lengths.size()) + " contigs, total " + number(total) + " bp, N50 " + number(n50) + " bp");
    });
  }
  void cleanup_directory(const std::string &path) {
    if (!directory(path)) return;
    for (const auto &name : entries(path)) {
      std::string entry = path + "/" + name;
      if (name == "flash" && opt.merge_pairs) continue;
      if (name == "flanking") { cleanup_directory(entry); if (entries(entry).empty()) remove_tree(entry); continue; }
      bool k_directory = name.size() > 1 && name[0] == 'k' &&
          name.find_first_not_of("0123456789", 1) == std::string::npos && directory(entry);
      if (k_directory) {
        for (const auto &child : entries(entry)) {
          bool retain = ends_with(child, ".contig_ends.fa") || ends_with(child, ".contig_ends.extended.fa") ||
              ends_with(child, ".contig_ends.extended.fa.info") || (opt.merge_pairs && child == "flank_mapping");
          if (!retain) remove_tree(entry + "/" + child);
        }
        if (entries(entry).empty()) remove_tree(entry);
      } else remove_tree(entry);
    }
  }
 public:
  Driver(Options options, std::string path) : opt(std::move(options)), executable(std::move(path)) {}
  ~Driver() {
    try { stop_producers(); if (!runtime_directory.empty()) remove_tree(runtime_directory); } catch (...) {}
  }
  int execute() {
    setup();
    if (opt.continue_mode && exists(opt.out_dir + "/done") && exists(final_path())) {
      log("Assembly already completed: " + final_path()); return 0;
    }
    merge_reads(); library();
    if (opt.legacy_k_list) log("Shared k list (legacy): " + Json(opt.k_list).dump());
    log("Reduction k list: " + Json(opt.reduction_k_list).dump());
    log("Flanking k list: " + Json(opt.flanking_k_list).dump());
    if (opt.flanking_consensus_only)
      log("Flanking consensus-only mode: retained contigs and consensus; no iterative read-derived edges");
    try {
      const auto &reduction = opt.reduction_k_list;
      int first = initial_k();
      bool first_built = opt.kmin_1pass && reduction.empty();
      if (first_built) {
        std::string output = graph(first);
        checkpoint("first-one-pass", {output + ".sdbg_info"}, [&]() {
          Args args{"read2sdbg", "-k", number(first), "-m", number(opt.min_count),
                    "--host_mem", number(host_mem), "--mem_flag", number(opt.mem_flag),
                    "--output_prefix", output, "--num_cpu_threads", number(opt.num_cpu_threads), "--read_lib_file", lib()};
          if (!opt.no_mercy) args.push_back("--need_mercy");
          core(args, "Build first graph in one pass");
        });
      } else count(first, true);
      for (size_t i = 1; i < reduction.size(); ++i) {
        int from = reduction[i - 1], to = reduction[i];
        count(to, false);
        // Resolve the input only when this stage runs: on resume a previously
        // generated reduced prefix must not accidentally become its own input.
        std::string reduced;
        if (opt.continue_mode && completed.value("reduce-k" + number(from) + "-k" + number(to), false)
            && resume_valid && exists(reduction_path(from, "reduced_edges") + ".edges.info")) {
          reduced = reduction_path(from, "reduced_edges");
          visited["reduce-k" + number(from) + "-k" + number(to)] = true;
          log("Reuse reduction edges for k = " + number(from));
        } else reduced = reduce(from, to);
        std::string assembled = reduction_path(from, "assembly_graph");
        graph_edges(from, reduced, assembled, "reduced-graph-k" + number(from));
        assemble(from, assembled, to); merge_edges(from, to);
      }
      int current;
      if (!reduction.empty()) {
        current = reduction.back();
        if (reduction.size() == 1) { graph_sequences(current, 0); assemble(current, graph(current)); }
        else {
          std::string output = reduction_path(current, "assembly_graph");
          graph_edges(current, current_edges(current), output, "handoff-graph-k" + number(current));
          assemble(current, output);
        }
        if (!opt.no_local) local(current, opt.flanking_k_list.front());
        std::string source = contig(current); int source_k = current;
        flanking = true; mkdirs(contig_dir());
        log("Restarting consensus flanking assembly at k = " + number(opt.flanking_k_list.front()) + " from k = " + number(source_k));
        seed(source_k, source); current = opt.flanking_k_list.front(); assemble(current, graph(current));
      } else {
        current = first;
        if (!first_built) graph_sequences(current, 0);
        assemble(current, graph(current));
      }
      for (int next : opt.flanking_k_list) if (next > current) {
        if (!opt.no_local) local(current, next);
        if (!opt.flanking_consensus_only) iterate(current, next);
        graph_sequences(next, current); assemble(next, graph(next)); current = next;
      }
      merge_final(current);
    } catch (const EarlyTerminate &early) { merge_final(early.k); }
    if (!opt.keep_tmp_files) {
      if (absolute(opt.temp_dir) == absolute(opt.out_dir + "/tmp")) cleanup_directory(opt.temp_dir);
      else remove_tree(opt.temp_dir);
    }
    atomic_text(opt.out_dir + "/done", "");
    log("ALL DONE. Time elapsed: " + number(std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count()) + " seconds");
    if (opt.test_mode && !opt.keep_tmp_files) remove_tree(opt.out_dir);
    return 0;
  }
};
}  // namespace

int workflow_main(int argc, char **argv) {
  try {
    install_signals();
    std::string self = executable_path(argv[0]);
    Options options = parse(argc, argv);
    if (options.test_mode && options.out_dir.empty()) {
      const char *temp = std::getenv("TMPDIR");
      options.out_dir = temporary_directory(std::string(temp ? temp : "/tmp") + "/megahit-test-");
      options.force_overwrite = true;
    }
    Driver driver(std::move(options), self);
    return driver.execute();
  } catch (const std::exception &error) { std::cerr << "megahit: " << error.what() << '\n'; return 1; }
}
}  // namespace native
