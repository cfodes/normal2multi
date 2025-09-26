#pragma once
#include <vector>
#include <string>
#include <stdexcept>
#include <Eigen/Core>



// ========================
// 点类定义
// ========================
template<typename _Ty>
//class Point : public std::vector<_Ty>
class Point
{
public:
    // 默认构造函数, 将所有元素设置为dat
    explicit Point(_Ty const& xx = _Ty(), _Ty const& yy = _Ty(), _Ty const& zz = _Ty())
    {
        x = xx;
        y = yy;
        z = zz;
    }
    typedef _Ty value_type;
    typedef unsigned short int size_type;
    _Ty x;
    _Ty y;
    _Ty z;
    // 复制构造函数
    inline Point(Point<_Ty> const& other)
    {
        (*this).x = other.x;
        (*this).y = other.y;
        (*this).z = other.z;
    }
    inline _Ty& operator[] (int const& i)
    {
        if (i == 0)
        {
            return (*this).x;
        }
        if (i == 1)
        {
            return (*this).y;
        }
        if (i == 2)
        {
            return (*this).z;
        }
    }


    inline _Ty const& operator[] (int const& i) const
    {

        if (i == 0)
        {
            return (*this).x;
        }
        if (i == 1)
        {
            return (*this).y;
        }
        if (i == 2)
        {
            return (*this).z;
        }
    }

    // operator =
    inline Point<_Ty> operator = (Point<_Ty> const& R)
    {
        (*this).x = R.x;
        (*this).y = R.y;
        (*this).z = R.z;

        return *this;
    }
    // operator -=
    inline void operator -= (Point<_Ty> const& R)
    {
        (*this).x -= R.x;
        (*this).y -= R.y;
        (*this).z -= R.z;
    }
    inline int const size() const
    {
        return 3;
    }
    inline size_type size()
    {
        return 3;
    }
    // 点之平方
    inline _Ty square() const
    {
        _Ty Re = _Ty();
        Re += (*this).x * (*this).x;
        Re += (*this).y * (*this).y;
        Re += (*this).z * (*this).z;

        return Re;
    }
    inline friend _Ty square(Point<_Ty> const& p)
    {
        _Ty Re = _Ty();

        Re += p.x * p.x;
        Re += p.y * p.y;
        Re += p.z * p.z;

        return Re;
    }
};


// ========================
// 节点类
// ========================
class Node
{
public:
    Point<double> point;
    Eigen::Vector3d df;
    int id;

    // 构造函数，所有参数都有默认值
    explicit Node(Point<double> p = Point<double>(), int i_d = 0, Eigen::Vector3d d_f = Eigen::Vector3d::Zero())
        : point(p), id(i_d), df(d_f) {}

    // 拷贝构造函数
    inline Node(Node const& other) : point(other.point), id(other.id), df(other.df) {}
};

// ========================
// 网格单元类
// ========================
//网格单元类定义
class element
{
public:
    int e_type;  // 单元类型
    int e_id;    // 单元编号
    std::vector<int> node_id;  // 节点编号数组（动态大小）

    // 构造函数，根据 e_type 初始化
    explicit element(int e_type_input)
    {
        e_type = e_type_input;
        e_id = 0;

        // 根据单元类型初始化节点数
        switch (e_type)
        {
        case 3:  // Line
            node_id.resize(2);  // 线单元有2个节点
            break;
        case 5:  // Triangle
            node_id.resize(3);  // 三角形单元有3个节点
            break;
        case 9:  // Quad
            node_id.resize(4);  // 四边形单元有4个节点
            break;
        case 10: // Tetrahedral
            node_id.resize(4);  // 四面体单元有4个节点
            break;
        case 12: // Hexahedral
            node_id.resize(8);  // 六面体单元有8个节点
            break;
        case 13: // Prism
            node_id.resize(6);  // 棱锥单元有6个节点
            break;
        case 14: // Pyramid
            node_id.resize(5);  // 金字塔单元有5个节点
            break;
        default:
            throw std::invalid_argument("Invalid element type");  // 如果e_type不合法，抛出异常
        }
    }
};

// ========================
// 边界类
// ========================
class boundary
{
public:
    int num_bound_elems;   //该边界单元数
    std::string boundary_tag;   //边界命名
    std::vector<element> bound_elements;   //边界单元构成的vector数组
};