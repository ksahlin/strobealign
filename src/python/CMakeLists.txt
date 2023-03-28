find_package(Python 3.10 COMPONENTS Interpreter Development.Module REQUIRED)
execute_process(
COMMAND "${Python_EXECUTABLE}" -c "import nanobind; print(nanobind.cmake_dir())"
OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NB_DIR)
list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(strobealign_python strobealign.cpp)
target_link_libraries(strobealign_python PRIVATE salib)
set_property(TARGET strobealign_python PROPERTY OUTPUT_NAME strobealign)
