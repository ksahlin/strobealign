#ifndef TMPDIR_HPP
#define TMPDIR_HPP

#include <filesystem>
#include <sstream>
#include <string>

/* Temporary directory that gets deleted automatically by the destructor */
class TemporaryDirectory {
public:
    TemporaryDirectory(const std::string& prefix = "tmp") {
        int i = 1;
        std::filesystem::path dir;
        while (true) {
            std::stringstream s;
            s << prefix << "-" << i;
            dir = std::filesystem::temp_directory_path() / s.str();
            if (std::filesystem::create_directory(dir)) {
                break;
            }
            i++;
        }
        m_path = dir;
    }

    ~TemporaryDirectory() {
        std::filesystem::remove_all(m_path);
    }

    const std::filesystem::path& path() const {
        return m_path;
    }

private:
    std::filesystem::path m_path;
};

#endif
