#ifndef MEGAHIT_NATIVE_PLATFORM_H
#define MEGAHIT_NATIVE_PLATFORM_H

#include <cstdint>
#include <functional>
#include <iosfwd>
#include <stdexcept>
#include <string>
#include <vector>
#include "vendor/nlohmann/json.hpp"

namespace native {
using Json = nlohmann::json;
struct Error : std::runtime_error { using std::runtime_error::runtime_error; };
bool exists(const std::string &path);
bool directory(const std::string &path);
uint64_t size(const std::string &path);
std::string parent(const std::string &path);
std::string base(const std::string &path);
std::string absolute(const std::string &path);
bool ends_with(const std::string &value, const std::string &suffix);
void mkdirs(const std::string &path);
void remove_file(const std::string &path);
void remove_tree(const std::string &path);
std::vector<std::string> entries(const std::string &path);
void replace(const std::string &source, const std::string &dest);
void copy(const std::string &source, const std::string &dest);
void atomic_text(const std::string &path, const std::string &text);
Json read_json(const std::string &path);
void atomic_json(const std::string &path, const Json &value);
Json signature(const std::string &path);
std::string executable_path(const char *argv0);
std::string temporary_directory(const std::string &prefix);
void decompress(const std::string &source, const std::string &output);
void convert_fastq(const std::string &source, const std::string &output,
                   const std::string &pigz, int threads,
                   const std::function<void(const std::string &)> &log);
void read_fasta(const std::string &path,
                const std::function<void(const std::string &, const std::string &)> &record);
void run(const std::string &executable, const std::vector<std::string> &arguments,
         const std::function<void(const std::string &)> &stderr_line,
         const std::string &stdout_path = "");
void install_signals();
extern const unsigned char bundled_pigz[];
extern const size_t bundled_pigz_size;
}  // namespace native
#endif
