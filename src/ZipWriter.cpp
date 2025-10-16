#include "ZipWriter.hpp"

#include <array>
#include <fstream>
#include <stdexcept>

namespace {
constexpr uint32_t CRC32_POLY = 0xEDB88320u;

std::array<uint32_t, 256> make_crc_table()
{
    std::array<uint32_t, 256> table{};
    for (uint32_t i = 0; i < 256; ++i) {
        uint32_t c = i;
        for (int j = 0; j < 8; ++j) {
            if (c & 1) {
                c = CRC32_POLY ^ (c >> 1);
            } else {
                c >>= 1;
            }
        }
        table[i] = c;
    }
    return table;
}

const std::array<uint32_t, 256>& crc_table()
{
    static const std::array<uint32_t, 256> table = make_crc_table();
    return table;
}
} // namespace

void ZipWriter::add_file(const std::string& path, const std::string& contents)
{
    entries_.push_back(Entry{path, contents});
}

void ZipWriter::save(const std::string& filename) const
{
    std::ofstream ofs(filename, std::ios::binary | std::ios::trunc);
    if (!ofs) {
        throw std::runtime_error("Failed to open zip file: " + filename);
    }

    struct CentralDirInfo {
        std::string path;
        uint32_t crc;
        uint32_t compressed_size;
        uint32_t uncompressed_size;
        uint32_t local_header_offset;
    };

    std::vector<CentralDirInfo> central_directory;
    central_directory.reserve(entries_.size());

    uint32_t offset = 0;

    for (const auto& entry : entries_) {
        const uint32_t crc = crc32(entry.contents);
        const uint32_t size = static_cast<uint32_t>(entry.contents.size());
        const uint16_t name_len = static_cast<uint16_t>(entry.path.size());

        // Local file header
        write_le32(ofs, 0x04034B50u); // signature
        write_le16(ofs, 20);          // version needed
        write_le16(ofs, 0);           // general purpose bit flag
        write_le16(ofs, 0);           // compression method: stored
        write_le16(ofs, 0);           // last mod time
        write_le16(ofs, 0);           // last mod date
        write_le32(ofs, crc);
        write_le32(ofs, size);
        write_le32(ofs, size);
        write_le16(ofs, name_len);
        write_le16(ofs, 0);           // extra length
        ofs.write(entry.path.data(), name_len);
        ofs.write(entry.contents.data(), entry.contents.size());

        CentralDirInfo info;
        info.path = entry.path;
        info.crc = crc;
        info.compressed_size = size;
        info.uncompressed_size = size;
        info.local_header_offset = offset;
        central_directory.push_back(info);

        offset += 30 + name_len + size; // local header fixed 30 bytes + name + data
    }

    const uint32_t central_dir_offset = offset;

    for (const auto& info : central_directory) {
        const uint16_t name_len = static_cast<uint16_t>(info.path.size());

        write_le32(ofs, 0x02014B50u); // central directory signature
        write_le16(ofs, 20);          // version made by
        write_le16(ofs, 20);          // version needed
        write_le16(ofs, 0);           // general purpose bit flag
        write_le16(ofs, 0);           // compression method
        write_le16(ofs, 0);           // mod time
        write_le16(ofs, 0);           // mod date
        write_le32(ofs, info.crc);
        write_le32(ofs, info.compressed_size);
        write_le32(ofs, info.uncompressed_size);
        write_le16(ofs, name_len);
        write_le16(ofs, 0);           // extra length
        write_le16(ofs, 0);           // file comment length
        write_le16(ofs, 0);           // disk number start
        write_le16(ofs, 0);           // internal file attributes
        write_le32(ofs, 0);           // external file attributes
        write_le32(ofs, info.local_header_offset);
        ofs.write(info.path.data(), name_len);

        offset += 46 + name_len; // central directory entry size
    }

    const uint32_t central_dir_size = offset - central_dir_offset;

    write_le32(ofs, 0x06054B50u); // end of central dir signature
    write_le16(ofs, 0);           // number of this disk
    write_le16(ofs, 0);           // number of the disk with the start of the central directory
    write_le16(ofs, static_cast<uint16_t>(central_directory.size()));
    write_le16(ofs, static_cast<uint16_t>(central_directory.size()));
    write_le32(ofs, central_dir_size);
    write_le32(ofs, central_dir_offset);
    write_le16(ofs, 0);           // comment length

    ofs.close();
}

uint32_t ZipWriter::crc32(const std::string& data)
{
    uint32_t c = 0xFFFFFFFFu;
    for (unsigned char ch : data) {
        c = crc_table()[(c ^ ch) & 0xFFu] ^ (c >> 8);
    }
    return c ^ 0xFFFFFFFFu;
}

void ZipWriter::write_le16(std::ostream& os, uint16_t value)
{
    char bytes[2];
    bytes[0] = static_cast<char>(value & 0xFFu);
    bytes[1] = static_cast<char>((value >> 8) & 0xFFu);
    os.write(bytes, 2);
}

void ZipWriter::write_le32(std::ostream& os, uint32_t value)
{
    char bytes[4];
    bytes[0] = static_cast<char>(value & 0xFFu);
    bytes[1] = static_cast<char>((value >> 8) & 0xFFu);
    bytes[2] = static_cast<char>((value >> 16) & 0xFFu);
    bytes[3] = static_cast<char>((value >> 24) & 0xFFu);
    os.write(bytes, 4);
}
