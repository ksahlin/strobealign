#include "version.hpp"

std::string version_string() {
    if (git_IsPopulated()) {
        return std::string(git_Describe());
    } else {
        return VERSION_STRING;
    }
}
