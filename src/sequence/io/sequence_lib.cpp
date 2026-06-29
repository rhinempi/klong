//
// Created by vout on 6/29/19.
//

#include "sequence_lib.h"
#include "async_sequence_reader.h"

#include <algorithm>
#include <cstdio>
#include <deque>
#include <future>
#include <stdexcept>

namespace {

struct LibSpec {
  size_t id;
  std::string metadata;
  std::string type;
  std::string file_name1;
  std::string file_name2;
  std::string temp_file;
};

struct LibBuildResult {
  size_t id;
  std::string metadata;
  std::string type;
  std::string temp_file;
  int64_t reads;
  int64_t bases;
  unsigned max_read_len;
};

std::unique_ptr<BaseSequenceReader> MakeReader(const LibSpec &spec) {
  if (spec.type == "pe") {
    return std::unique_ptr<BaseSequenceReader>(
        new PairedFastxReader(spec.file_name1, spec.file_name2));
  }
  if (spec.type == "se" || spec.type == "interleaved") {
    return std::unique_ptr<BaseSequenceReader>(new FastxReader(spec.file_name1));
  }
  throw std::runtime_error("Cannot identify read library type " + spec.type);
}

LibBuildResult BuildOneLib(const LibSpec &spec, bool overlap_read_write) {
  std::ofstream bin_file(spec.temp_file,
                         std::ofstream::binary | std::ofstream::out);
  if (!bin_file.is_open()) {
    throw std::runtime_error("Cannot open temporary read library " +
                             spec.temp_file);
  }

  auto reader = MakeReader(spec);

  int64_t reads = 0;
  int64_t bases = 0;
  unsigned max_read_len = 0;

  if (overlap_read_write) {
    AsyncSequenceReader async_reader(reader.get());

    while (true) {
      const auto &seq_batch = async_reader.Next();
      if (seq_batch.seq_count() == 0) {
        break;
      }

      reads += seq_batch.seq_count();
      bases += seq_batch.base_count();
      seq_batch.WriteSequences(bin_file);
      max_read_len = std::max(max_read_len, seq_batch.max_length());
    }
  } else {
    SeqPackage seq_batch;

    while (true) {
      seq_batch.Clear();
      reader->Read(&seq_batch, AsyncSequenceReader::kDefaultNumSeqPerBatch,
                   AsyncSequenceReader::kDefaultNumBasesPerBatch, false);
      if (seq_batch.seq_count() == 0) {
        break;
      }

      reads += seq_batch.seq_count();
      bases += seq_batch.base_count();
      seq_batch.WriteSequences(bin_file);
      max_read_len = std::max(max_read_len, seq_batch.max_length());
    }
  }

  if (!bin_file.good()) {
    throw std::runtime_error("Failed while writing temporary read library " +
                             spec.temp_file);
  }

  if ((spec.type == "pe" || spec.type == "interleaved") && reads % 2 != 0) {
    throw std::runtime_error("PE library number of reads is odd: " +
                             std::to_string(reads) + " for " +
                             spec.metadata);
  }

  return {spec.id, spec.metadata, spec.type, spec.temp_file, reads, bases,
          max_read_len};
}

}  // namespace

void SequenceLibCollection::Build(const std::string &lib_file,
                                  const std::string &out_prefix,
                                  unsigned num_threads) {
  std::ifstream lib_config(lib_file);

  if (!lib_config.is_open()) {
    xfatal("File to open read_lib file: {}\n", lib_file.c_str());
  }

  std::string metadata;
  std::string type;
  std::string file_name1;
  std::string file_name2;

  std::vector<LibSpec> specs;

  while (std::getline(lib_config, metadata)) {
    if (!(lib_config >> type)) {
      xfatal("Cannot read read library type for {s}\n", metadata.c_str());
    }

    file_name1.clear();
    file_name2.clear();
    if (type == "pe") {
      lib_config >> file_name1 >> file_name2;
    } else if (type == "se" || type == "interleaved") {
      lib_config >> file_name1;
    } else {
      xerr("Cannot identify read library type {}\n", type.c_str());
      xfatal("Valid types: pe, se, interleaved\n");
    }

    specs.push_back({specs.size(), metadata, type, file_name1, file_name2,
                     out_prefix + ".part." + std::to_string(specs.size()) +
                         ".bin"});
    std::getline(lib_config, metadata);  // eliminate the "\n"
  }

  if (specs.empty()) {
    xfatal("No read library found in {s}\n", lib_file.c_str());
  }

  unsigned worker_count = std::max(1u, num_threads);
  worker_count = std::min(worker_count, static_cast<unsigned>(specs.size()));

  xinfo("Converting {} read libraries using {} worker(s)\n",
        specs.size(), worker_count);

  int64_t total_reads = 0;
  int64_t total_bases = 0;

  std::vector<SequenceLib> libs;

  std::ofstream bin_file(out_prefix + ".bin",
                         std::ofstream::binary | std::ofstream::out);
  if (!bin_file.is_open()) {
    xfatal("Cannot open binary read library {s}\n",
           (out_prefix + ".bin").c_str());
  }

  std::deque<std::future<LibBuildResult>> futures;
  size_t next_spec = 0;

  auto launch_next = [&]() {
    const auto spec = specs[next_spec++];
    const bool overlap_read_write = worker_count == 1;
    futures.emplace_back(std::async(std::launch::async,
                                    [spec, overlap_read_write]() {
                                      return BuildOneLib(spec,
                                                         overlap_read_write);
                                    }));
  };

  auto append_result = [&](const LibBuildResult &result) {
    int64_t begin_index = total_reads;
    total_reads += result.reads;
    total_bases += result.bases;

    std::ifstream part_file(result.temp_file, std::ifstream::binary);
    if (!part_file.is_open()) {
      xfatal("Cannot open temporary read library {s}\n",
             result.temp_file.c_str());
    }
    bin_file << part_file.rdbuf();
    if (!bin_file.good()) {
      xfatal("Failed while writing binary read library {s}\n",
             (out_prefix + ".bin").c_str());
    }
    part_file.close();
    std::remove(result.temp_file.c_str());

    xinfo("Lib {} ({s}): {s}, {} reads, {} max length\n",
          result.id, result.metadata.c_str(), result.type.c_str(),
          result.reads, result.max_read_len);

    libs.emplace_back(nullptr, begin_index, total_reads, result.max_read_len,
                      result.type != "se", result.metadata);
  };

  try {
    while (next_spec < specs.size() || !futures.empty()) {
      while (next_spec < specs.size() && futures.size() < worker_count) {
        launch_next();
      }

      auto result = futures.front().get();
      futures.pop_front();
      append_result(result);
    }
  } catch (const std::exception &e) {
    xfatal("Failed to build read library: {s}\n", e.what());
  }
  bin_file.close();

  std::ofstream lib_info_file(out_prefix + ".lib_info");
  lib_info_file << total_bases << ' ' << total_reads << '\n';

  for (auto &lib : libs) {
    lib.DumpMetadata(lib_info_file);
  }
  lib_info_file.close();
}

void SequenceLibCollection::Read(SeqPackage *pkg, bool reverse_seq) {
  std::ifstream lib_info_file(path_ + ".lib_info");
  int64_t total_bases, num_reads;
  bool is_paired;
  std::string metadata;

  lib_info_file >> total_bases >> num_reads;
  std::getline(lib_info_file, metadata);  // eliminate the "\n"

  while (std::getline(lib_info_file, metadata)) {
    int64_t start, end;
    int max_read_len;
    lib_info_file >> start >> end >> max_read_len >> is_paired;
    libs_.emplace_back(pkg, start, end, max_read_len, is_paired, metadata);
    std::getline(lib_info_file, metadata);  // eliminate the "\n"
  }

  pkg->Clear();
  pkg->ReserveSequences(num_reads);
  pkg->ReserveBases(total_bases);
  BinaryReader reader(path_ + ".bin");

  xinfo("Before reading, sizeof seq_package: {}\n", pkg->size_in_byte());
  reader.ReadAll(pkg, reverse_seq);
  xinfo("After reading, sizeof seq_package: {}\n", pkg->size_in_byte());
}

std::pair<int64_t, int64_t> SequenceLibCollection::GetSize() const {
  std::ifstream lib_info_file(path_ + ".lib_info");
  int64_t total_bases, num_reads;
  lib_info_file >> total_bases >> num_reads;
  return {total_bases, num_reads};
}
