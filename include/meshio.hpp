#pragma once

#include <string>
#include "State.hpp"

// 读取网格文件
void readfile(std::string const& file_name, State& S);

// 写网格文件
void writefile(std::string const& file_name, const State& S);

// 将物面节点写出
void write_wall(const std::vector<Node>& wall_nodes, const std::string& file_name = "wall.dat");

//字符串切割函数
//该函数实现将字符串str里按split分隔符分开的字符，分别存入res里
inline void Stringsplit(std::string str, const char split, std::vector<std::string>& res)
{
    std::istringstream iss(str); //将str改为输入流
    std::string token;           //接收缓冲区
    while (std::getline(iss, token, split))
    {
        res.push_back(token);
    }
}
//该函数实现将str按spilt分隔符进行分割，同时它不会存入空字符串并且适用于较大的str
inline std::vector<std::string> stringSplit(const std::string& str, char spilt)
{
    std::stringstream ss(str);
    std::string token;
    std::vector<std::string> elems;
    while (getline(ss, token, spilt))
    {
        if (!token.empty())
        {
            elems.push_back(token);
        }
    }
    return elems;
}