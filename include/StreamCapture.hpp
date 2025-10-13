#pragma once

#include <ostream>
#include <sstream>
#include <string>

/**
 * StreamCapture temporarily redirects an output stream (e.g. std::cout)
 * to an internal buffer. Once the work is done, the buffered content can be
 * restored to the original stream and/or flushed to a file without impacting
 * timing in the hot loop.
 */
class StreamCapture {
public:
    explicit StreamCapture(std::ostream& stream);
    ~StreamCapture();

    // Restore the original stream buffer immediately. Safe to call multiple times.
    void stop();

    // Return buffered output as a std::string.
    std::string str() const;

    // Append buffered output to another stream.
    void dump_to_stream(std::ostream& out) const;

    // Save buffered output to a file. Throws std::runtime_error on failure.
    void save_to_file(const std::string& file_path) const;

private:
    std::ostream& stream_;
    std::streambuf* original_buf_;
    std::ostringstream buffer_;
    bool active_;
};
