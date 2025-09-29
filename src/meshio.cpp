#include "meshio.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

//网格文件读取函数
void readfile(std::string const& file_name, State& S)
{
    std::ifstream ifs;
    ifs.open(file_name, std::ios::in);

    if (!ifs.is_open())
    {
        std::cout << "文件打开失败！" << std::endl;
        return;
    }
    int Ndime = 0;  //网格维度
    int Npoint = 0; //节点数目
    int Nelem = 0;  //单元数目
    int NBOUND = 0; //边界的数目
    int N_subbound_elems = 0;  //子边界的单元数
    int i_boundary = 0; //用于确定读到第几个边界，i_boundary=1代表读到第一个单元
    std::string BMARKER;  //边界的命名

    int ielem_type = 0; //读取到的单元类型
    int line = 0;
    bool judge_start_read_elem = false;
    bool judge_start_read_point = false;
    bool judge_start_read_Belem = false;

    std::string buf;     //每一行的字符串
    std::cout.precision(13);  //设置输出精度

    int step = 0;

    while (std::getline(ifs, buf))
    {
        step += 1;
        if (buf.find("%") != std::string::npos) //一旦读到%，就将所有读取信息的功能关闭
        {
            judge_start_read_elem = false;
            judge_start_read_point = false;
            judge_start_read_Belem = false;
        }
        if (buf.find("NDIME") != std::string::npos)  //找到NDIME= 2(3)，将2提取出来
        {
            std::vector<std::string> s_NDIME = stringSplit(buf, ' ');
            if (s_NDIME.size() > 1)
            {
                std::stringstream(s_NDIME[1]) >> Ndime;
            }
            else//防止NDIME= 2(3)中间没有空格
            {
                std::vector<std::string> s_NDIME1 = stringSplit(buf, '=');
                std::stringstream(s_NDIME1[1]) >> Ndime;
            }
            S.Dimension = Ndime;
            
            std::cout << "=====================" << std::endl;
            std::cout << " Mesh Information " << std::endl;
            std::cout << "=====================" << std::endl;
            std::cout << " Mesh Dimension     : " << S.Dimension << std::endl;
       
            
        }
        if (buf.find("NELEM") != std::string::npos)
        {
            std::vector<std::string> s_NELE = stringSplit(buf, ' ');
            if (s_NELE.size() > 1)
            {
                std::stringstream(s_NELE[1]) >> Nelem;
            }
            else//防止NDIME= 2(3)中间没有空格
            {
                std::vector<std::string> s_NELE1 = stringSplit(buf, '=');
                std::stringstream(s_NELE1[1]) >> Nelem;
            }
            std::cout << "The number of the elements is " << Nelem << std::endl;
            judge_start_read_elem = true; //开始读取单元信息
            S.Nelements = Nelem;

            std::cout << " Number of Elements : " << Nelem << std::endl;
            std::cout << " Status             : Reading elements..." << std::endl;

            continue;
        }
        if (buf.find("NPOIN") != std::string::npos)
        {
            std::vector<std::string> s_NPOIN = stringSplit(buf, ' ');
            if (s_NPOIN.size() > 1)
            {
                std::stringstream(s_NPOIN[1]) >> Npoint;
            }
            else//防止=后没有空格
            {
                std::vector<std::string> s_NPOIN1 = stringSplit(buf, '=');
                std::stringstream(s_NPOIN1[1]) >> Npoint;
            }
            S.NPoints = Npoint;

            std::cout << " Number of Points   : " << S.NPoints << std::endl;
            std::cout << " Status             : Reading points..." << std::endl;

            judge_start_read_point = true;   //开始读取节点信息
            
            continue;
        }
        if (buf.find("NMARK") != std::string::npos)
        {
            std::vector<std::string> s_NBOUND = stringSplit(buf, ' ');
            if (s_NBOUND.size() > 1)
            {
                std::stringstream(s_NBOUND[1]) >> NBOUND;
            }
            else//防止=后没有空格
            {
                std::vector<std::string> s_NBOUND1 = stringSplit(buf, '=');
                std::stringstream(s_NBOUND1[1]) >> NBOUND;
            }
            S.every_boundary.resize(NBOUND);

            std::cout << "=====================" << std::endl;
            std::cout << " Boundary Information" << std::endl;
            std::cout << "=====================" << std::endl;
            std::cout << " Number of Boundaries: " << NBOUND << std::endl;

            continue;
        }
        if (buf.find("MARKER_TAG") != std::string::npos)
        {
            ++i_boundary;  //读到TAG意味着读到一个新的边界
            std::vector<std::string> s_MAREKER = stringSplit(buf, ' ');
            if (s_MAREKER.size() > 1)
            {
                std::stringstream(s_MAREKER[1]) >> BMARKER;
            }
            else//防止=后没有空格
            {
                std::vector<std::string> s_MAREKER1 = stringSplit(buf, '=');
                std::stringstream(s_MAREKER1[1]) >> BMARKER;
            }
            S.every_boundary[i_boundary - 1].boundary_tag = BMARKER;
            judge_start_read_Belem = false;  //在这一行禁止开始读边界单元，因为下一行还是文本信息
            
            std::cout << " Boundary Tag       : " << S.every_boundary[i_boundary - 1].boundary_tag << std::endl;
            std::cout << " Status             : Reading boundary elements..." << std::endl;

            continue;
        }
        if (buf.find("MARKER_ELEMS") != std::string::npos)
        {
            std::vector<std::string> s_MAREKER_ELE = stringSplit(buf, ' ');
            if (s_MAREKER_ELE.size() > 1)
            {
                std::stringstream(s_MAREKER_ELE[1]) >> N_subbound_elems;
            }
            else//防止=后没有空格
            {
                std::vector<std::string> s_MAREKER_ELE1 = stringSplit(buf, '=');
                std::stringstream(s_MAREKER_ELE1[1]) >> N_subbound_elems;
            }
            S.every_boundary[i_boundary - 1].num_bound_elems = N_subbound_elems;
            judge_start_read_Belem = true;  //在这一行开始读边界单元，因为下一行开始是边界信息
            std::cout << "The number of the elements at this boundary is " << S.every_boundary[i_boundary - 1].num_bound_elems << std::endl;
            continue;
        }

        if (judge_start_read_elem == true)  //读完单元数后，开始读单元信息
        {
            //std::cout << step << endl;
            std::vector<std::string>s = stringSplit(buf, ' ');
            std::stringstream(s[0]) >> ielem_type;
            element e_i(ielem_type);
            for (int i = 0; i < s.size() - 2; ++i)
            {
                std::stringstream(s[i + 1]) >> e_i.node_id[i];
            }
            
            std::stringstream(s.back()) >> e_i.e_id;
            S.every_elements.push_back(e_i);
        }
        if (judge_start_read_point == true)  // 读完单元数后，开始读节点信息
        {
            std::vector<std::string> s = stringSplit(buf, ' ');
            if (s.size() < 3)
            {
                s = stringSplit(buf, '\t');
            }
            double pi_info[3]; // 三维就是三个坐标加一个 id，二维就是两个坐标加一个 id
            int pi_id = 0;

            // 读取 id 之前的坐标信息
            for (int i = 0; i < s.size() - 1; ++i)
            {
                std::stringstream(s[i]) >> pi_info[i];
            }
            // 读取 id
            std::stringstream(s.back()) >> pi_id;

            if (s.size() == 4) // 三维包括节点编号是3+1=4个数
            {
                Point<double> p_temp(pi_info[0], pi_info[1], pi_info[2]);
                Node pi(p_temp, pi_id);
                S.node_coords.push_back(pi);
            }
            else if (s.size() == 3) // 二维包括节点编号是2+1=3个数
            {
                Point<double> p_temp(pi_info[0], pi_info[1], 0);
                Node pi(p_temp, pi_id);
                S.node_coords.push_back(pi);
            }
        }
        if (judge_start_read_Belem == true)  //读完边界单元数，开始读边界单元
        {
            std::vector<std::string>s = stringSplit(buf, ' ');
            std::stringstream(s[0]) >> ielem_type;
            element e_i(ielem_type);
            for (int i = 0; i < s.size() - 1; ++i)
            {
                std::stringstream(s[i + 1]) >> e_i.node_id[i];
            }
            S.every_boundary[i_boundary - 1].bound_elements.push_back(e_i);

            /*调试代码
            std::cout << _every_boundary[i_boundary - 1].bound_elements.back().e_type << "\t";
            for (int i = 0; i < s.size() - 1; ++i)
            {
                std::cout << _every_boundary[i_boundary - 1].bound_elements.back().node_id[i] << "\t";
            }
            std::cout << std::endl;*/
        }
    }
}


//写网格文件函数
void writefile(const std::string& file_name, const State& S)
    
{
    std::ofstream ofs;
    ofs.open(file_name, std::ios::out | std::ios::trunc);
    if (!ofs.is_open())
    {
        std::cout << "文件打开失败！" << std::endl;
        return;
    }
    else
    {
        
        std::cout << "=====================" << std::endl;
        std::cout << " Writing mesh file..." << std::endl;
        std::cout << "=====================" << std::endl;
        std::cout << " file name: " << file_name << std::endl;

        //写维度
        ofs << "%" << std::endl;
        ofs << "% Problem dimension" << std::endl;
        ofs << "NDIME= " << S.Dimension << std::endl;
        ofs << "%" << std::endl;
        ofs << "% Inner element connectivity" << std::endl;
        ofs << "%" << std::endl;
        //写单元数
        ofs << "NELEM= " << S.Nelements << std::endl;
        //写单元
        std::for_each(S.every_elements.begin(), S.every_elements.end(), [&ofs](const auto& elem)
            {
                ofs << elem.e_type;
                for (const auto& node : elem.node_id)
                {
                    ofs << " " << node;
                }
                ofs << " " << elem.e_id;
                ofs << std::endl;
            });
        
        ofs << "%" << std::endl;
        ofs << "% Node coordinates" << std::endl;
        ofs << "%" << std::endl;
        //写节点数
        ofs << "NPOIN= " << S.NPoints << std::endl;
        //写节点
        std::for_each(S.node_coords.begin(), S.node_coords.end(), [&ofs, &S](const auto& node)
            {
                if (S.Dimension == 2)
                {
                    ofs << std::fixed << std::setprecision(13) << node.point.x << " " << node.point.y << " ";
                    ofs << node.id;
                    ofs << std::endl;
                }
                if (S.Dimension == 3)
                {
                    ofs << std::fixed << std::setprecision(13) << node.point.x << " " << node.point.y << " " << node.point.z << " ";
                    ofs << node.id;
                    ofs << std::endl;
                }
            });
        ofs << "%" << std::endl;
        ofs << "% Boundary elements" << std::endl;
        ofs << "%" << std::endl;
        //写边界数
        ofs << "NMARK= " << S.every_boundary.size() << std::endl;
        //写边界
        std::for_each(S.every_boundary.begin(), S.every_boundary.end(), [&ofs](const auto& iboundary)
        {
            ofs << "MARKER_TAG= " << iboundary.boundary_tag << std::endl;
            ofs << "MARKER_ELEMS= " << iboundary.num_bound_elems << std::endl;
            for (const auto& elem : iboundary.bound_elements)
            {
                ofs << elem.e_type;
                for (const auto& node : elem.node_id)
                {
                    ofs << " " << node;
                }
                ofs << std::endl;
            }
        });

        std::cout << "successfully write mesh!" << std::endl;
        std::cout << "=====================" << std::endl;
    }
}

