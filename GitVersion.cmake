execute_process(
  COMMAND git describe --tags --dirty
  OUTPUT_VARIABLE GIT_VERSION
  RESULT_VARIABLE GIT_ERROR
  ERROR_QUIET
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

if (NOT GIT_ERROR)
  string(SUBSTRING "${GIT_VERSION}" 1 -1 GIT_VERSION)  # Remove "v" prefix
  set(PROJECT_VERSION "${GIT_VERSION}")
else()
  set(PROJECT_VERSION "${DEFAULT_VERSION}")
endif()
message("Setting version to ${PROJECT_VERSION}")

configure_file("${CONFIGIN}" "${CONFIGOUT}")
