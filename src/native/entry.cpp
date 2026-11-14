#include "platform.h"
#include <cstring>
#include <iostream>

int megahit_core_main(int argc, char **argv);
namespace native { int workflow_main(int argc, char **argv); }

int main(int argc, char **argv) {
  if (argc > 1 && std::strcmp(argv[1], "--core") == 0)
    return megahit_core_main(argc - 1, argv + 1);
  if (argc > 1 && std::strcmp(argv[1], "--native-decompress") == 0) {
    if (argc != 4) return 2;
    try { native::decompress(argv[2], argv[3]); return 0; }
    catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 1; }
  }
  if (argc > 1 && std::strcmp(argv[1], "--native-info") == 0) {
    std::cout << native::Json({{"workflow", "native-cpp"}, {"python_required", false},
                              {"core_linked", true}, {"bundled_pigz", native::bundled_pigz_size > 0}}).dump() << '\n';
    return 0;
  }
  if (argc > 1 && argv[1][0] != '-') return megahit_core_main(argc, argv);
  return native::workflow_main(argc, argv);
}
