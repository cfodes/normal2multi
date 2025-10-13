#include "XlsxWriter.hpp"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>

XlsxCell XlsxCell::String(const std::string& value)
{
    XlsxCell cell;
    cell.type = Type::String;
    cell.text = value;
    return cell;
}

XlsxCell XlsxCell::Number(double value)
{
    XlsxCell cell;
    cell.type = Type::Number;
    cell.number = value;
    return cell;
}

void XlsxWriter::add_sheet(const std::string& name,
                           const std::vector<std::vector<XlsxCell>>& rows)
{
    std::vector<std::string> existing;
    existing.reserve(sheets_.size());
    for (const auto& s : sheets_) {
        existing.push_back(s.name);
    }

    Sheet sheet;
    sheet.name = sanitize_sheet_name(name, existing);
    sheet.rows = rows;
    sheets_.push_back(std::move(sheet));
}

void XlsxWriter::save(const std::string& filename) const
{
    if (sheets_.empty()) {
        throw std::runtime_error("XlsxWriter: no sheets to write");
    }

    ZipWriter zip;

    // Content Types
    std::ostringstream ct;
    ct << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    ct << "<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\">\n";
    ct << "  <Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/>\n";
    ct << "  <Default Extension=\"xml\" ContentType=\"application/xml\"/>\n";
    ct << "  <Override PartName=\"/xl/workbook.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml\"/>\n";
    ct << "  <Override PartName=\"/docProps/core.xml\" ContentType=\"application/vnd.openxmlformats-package.core-properties+xml\"/>\n";
    ct << "  <Override PartName=\"/docProps/app.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.extended-properties+xml\"/>\n";
    for (size_t i = 0; i < sheets_.size(); ++i) {
        ct << "  <Override PartName=\"/xl/worksheets/sheet" << (i + 1)
           << ".xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>\n";
    }
    ct << "</Types>\n";
    zip.add_file("[Content_Types].xml", ct.str());

    // Relationships root
    std::ostringstream rels;
    rels << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    rels << "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">\n";
    rels << "  <Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument\" Target=\"xl/workbook.xml\"/>\n";
    rels << "  <Relationship Id=\"rId2\" Type=\"http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties\" Target=\"docProps/core.xml\"/>\n";
    rels << "  <Relationship Id=\"rId3\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties\" Target=\"docProps/app.xml\"/>\n";
    rels << "</Relationships>\n";
    zip.add_file("_rels/.rels", rels.str());

    // workbook rels
    std::ostringstream workbook_rels;
    workbook_rels << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    workbook_rels << "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">\n";
    for (size_t i = 0; i < sheets_.size(); ++i) {
        workbook_rels << "  <Relationship Id=\"rId" << (i + 1)
                      << "\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" Target=\"worksheets/sheet"
                      << (i + 1) << ".xml\"/>\n";
    }
    workbook_rels << "</Relationships>\n";
    zip.add_file("xl/_rels/workbook.xml.rels", workbook_rels.str());

    // workbook xml
    std::ostringstream workbook;
    workbook << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n";
    workbook << "<workbook xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" xmlns:r=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships\">\n";
    workbook << "  <sheets>\n";
    for (size_t i = 0; i < sheets_.size(); ++i) {
        workbook << "    <sheet name=\"" << escape_xml(sheets_[i].name)
                 << "\" sheetId=\"" << (i + 1)
                 << "\" r:id=\"rId" << (i + 1) << "\"/>\n";
    }
    workbook << "  </sheets>\n";
    workbook << "</workbook>\n";
    zip.add_file("xl/workbook.xml", workbook.str());

    // docProps/core.xml
    std::ostringstream core;
    core << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n";
    core << "<cp:coreProperties xmlns:cp=\"http://schemas.openxmlformats.org/package/2006/metadata/core-properties\" xmlns:dc=\"http://purl.org/dc/elements/1.1/\" xmlns:dcterms=\"http://purl.org/dc/terms/\" xmlns:dcmitype=\"http://purl.org/dc/dcmitype/\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";
    core << "  <dc:creator>normal2multi</dc:creator>\n";
    core << "  <cp:lastModifiedBy>normal2multi</cp:lastModifiedBy>\n";
    core << "</cp:coreProperties>\n";
    zip.add_file("docProps/core.xml", core.str());

    // docProps/app.xml
    std::ostringstream app;
    app << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n";
    app << "<Properties xmlns=\"http://schemas.openxmlformats.org/officeDocument/2006/extended-properties\" xmlns:vt=\"http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes\">\n";
    app << "  <Application>normal2multi</Application>\n";
    app << "  <HeadingPairs>\n";
    app << "    <vt:vector size=\"2\" baseType=\"variant\">\n";
    app << "      <vt:variant><vt:lpstr>Worksheets</vt:lpstr></vt:variant>\n";
    app << "      <vt:variant><vt:i4>" << sheets_.size() << "</vt:i4></vt:variant>\n";
    app << "    </vt:vector>\n";
    app << "  </HeadingPairs>\n";
    app << "  <TitlesOfParts>\n";
    app << "    <vt:vector size=\"" << sheets_.size() << "\" baseType=\"lpstr\">\n";
    for (const auto& s : sheets_) {
        app << "      <vt:lpstr>" << escape_xml(s.name) << "</vt:lpstr>\n";
    }
    app << "    </vt:vector>\n";
    app << "  </TitlesOfParts>\n";
    app << "</Properties>\n";
    zip.add_file("docProps/app.xml", app.str());

    // worksheets
    for (size_t i = 0; i < sheets_.size(); ++i) {
        zip.add_file("xl/worksheets/sheet" + std::to_string(i + 1) + ".xml",
                     build_sheet_xml(sheets_[i]));
    }

    zip.save(filename);
}

std::string XlsxWriter::sanitize_sheet_name(const std::string& name,
                                            const std::vector<std::string>& existing)
{
    std::string sanitized;
    sanitized.reserve(name.size());
    const std::string forbidden = "[]:*?/\\";

    for (char ch : name) {
        if (static_cast<unsigned char>(ch) < 0x20 || forbidden.find(ch) != std::string::npos) {
            sanitized.push_back('_');
        } else {
            sanitized.push_back(ch);
        }
    }
    if (sanitized.empty()) {
        sanitized = "Sheet";
    }
    if (sanitized.size() > 31) {
        sanitized.resize(31);
    }

    std::string base = sanitized;
    int suffix = 1;
    auto exists = [&](const std::string& candidate) {
        return std::find(existing.begin(), existing.end(), candidate) != existing.end();
    };

    while (exists(sanitized)) {
        std::string extra = "_" + std::to_string(suffix++);
        size_t max_len = 31;
        if (extra.size() >= max_len) {
            sanitized = extra.substr(extra.size() - max_len);
        } else {
            sanitized = base;
            if (sanitized.size() + extra.size() > max_len) {
                sanitized.resize(max_len - extra.size());
            }
            sanitized += extra;
        }
    }

    return sanitized;
}

std::string XlsxWriter::escape_xml(const std::string& text)
{
    std::string result;
    result.reserve(text.size());
    for (char ch : text) {
        switch (ch) {
        case '&': result += "&amp;"; break;
        case '<': result += "&lt;"; break;
        case '>': result += "&gt;"; break;
        case '\"': result += "&quot;"; break;
        case '\'': result += "&apos;"; break;
        default: result.push_back(ch); break;
        }
    }
    return result;
}

std::string XlsxWriter::column_name(size_t index)
{
    std::string name;
    size_t n = index;
    do {
        char ch = static_cast<char>('A' + (n % 26));
        name.insert(name.begin(), ch);
        n = n / 26;
        if (n == 0) {
            break;
        }
        --n;
    } while (true);
    return name;
}

std::string XlsxWriter::build_sheet_xml(const Sheet& sheet)
{
    std::ostringstream ss;
    ss << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n";
    ss << "<worksheet xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\">\n";
    ss << "  <sheetData>\n";

    for (size_t row_idx = 0; row_idx < sheet.rows.size(); ++row_idx) {
        const auto& row = sheet.rows[row_idx];
        if (row.empty()) {
            continue;
        }
        ss << "    <row r=\"" << (row_idx + 1) << "\">";
        for (size_t col_idx = 0; col_idx < row.size(); ++col_idx) {
            const auto& cell = row[col_idx];
            const std::string cell_ref = column_name(col_idx) + std::to_string(row_idx + 1);
            if (cell.type == XlsxCell::Type::String) {
                ss << "<c r=\"" << cell_ref << "\" t=\"inlineStr\"><is><t>"
                   << escape_xml(cell.text) << "</t></is></c>";
            } else {
                std::ostringstream value;
                value.precision(15);
                value << cell.number;
                ss << "<c r=\"" << cell_ref << "\"><v>" << value.str() << "</v></c>";
            }
        }
        ss << "</row>\n";
    }

    ss << "  </sheetData>\n";
    ss << "</worksheet>\n";
    return ss.str();
}
