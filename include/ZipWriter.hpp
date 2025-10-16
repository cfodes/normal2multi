#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <ostream>

class ZipWriter {
public:
    void add_file(const std::string& path, const std::string& contents);
    void save(const std::string& filename) const;

private:
    struct Entry {
        std::string path;
        std::string contents;
    };

    static uint32_t crc32(const std::string& data);
    static void write_le16(std::ostream& os, uint16_t value);
    static void write_le32(std::ostream& os, uint32_t value);

    std::vector<Entry> entries_;
};
