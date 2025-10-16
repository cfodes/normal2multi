#include "StreamCapture.hpp"

#include <fstream>
#include <stdexcept>

StreamCapture::StreamCapture(std::ostream& stream)
    : stream_(stream),
      original_buf_(stream.rdbuf()),
      buffer_(),
      active_(true)
{
    stream_.flush();
    stream_.rdbuf(buffer_.rdbuf());
}

StreamCapture::~StreamCapture()
{
    stop();
}

void StreamCapture::stop()
{
    if (active_) {
        stream_.flush();
        stream_.rdbuf(original_buf_);
        active_ = false;
    }
}

std::string StreamCapture::str() const
{
    return buffer_.str();
}

void StreamCapture::dump_to_stream(std::ostream& out) const
{
    out << buffer_.str();
}

void StreamCapture::save_to_file(const std::string& file_path) const
{
    std::ofstream ofs(file_path, std::ios::out | std::ios::trunc);
    if (!ofs) {
        throw std::runtime_error("Failed to open log file: " + file_path);
    }
    ofs << buffer_.str();
}
