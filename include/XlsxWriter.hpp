#pragma once

#include "ZipWriter.hpp"
#include <string>
#include <vector>

struct XlsxCell {
    enum class Type { String, Number };

    static XlsxCell String(const std::string& value);
    static XlsxCell Number(double value);

    Type type;
    std::string text;
    double number = 0.0;
};

class XlsxWriter {
public:
    void add_sheet(const std::string& name,
                   const std::vector<std::vector<XlsxCell>>& rows);

    void save(const std::string& filename) const;

private:
    struct Sheet {
        std::string name;
        std::vector<std::vector<XlsxCell>> rows;
    };

    static std::string sanitize_sheet_name(const std::string& name,
                                           const std::vector<std::string>& existing);
    static std::string escape_xml(const std::string& text);
    static std::string column_name(size_t index);
    static std::string build_sheet_xml(const Sheet& sheet);

    std::vector<Sheet> sheets_;
};
