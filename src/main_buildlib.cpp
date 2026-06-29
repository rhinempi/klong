#include "sequence/io/sequence_lib.h"
#include "utils/utils.h"

#include <algorithm>

void DisplayHelp(const char *program) {
  pfprintf(stderr, "Usage {s} <read_lib_file> <out_prefix> [num_threads]\n",
           program);
}

int main_build_lib(int argc, char **argv) {
  AutoMaxRssRecorder recorder;

  if (argc < 3) {
    DisplayHelp(argv[0]);
    exit(1);
  }
  unsigned num_threads = 1;
  if (argc >= 4) {
    num_threads = std::max(1, atoi(argv[3]));
  }
  SequenceLibCollection::Build(argv[1], argv[2], num_threads);

  return 0;
}
