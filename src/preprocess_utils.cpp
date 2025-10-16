#include "preprocess_utils.hpp"

void _set_wall_tree(GridBTree<int, Point<double>>& global_wall_tree,
                    const std::vector<Node>& wall)
//构建物面节点的二叉树，用以查询空间节点到物面节点的最小距离
//输入：global_wall_tree：全局的物面二叉树
//      wall：全局的物面节点数组
{
    global_wall_tree.setGridSize(static_cast<int>(wall.size()));    
    for (int i = 0; i < static_cast<int>(wall.size()); ++i) {
        global_wall_tree.setGrid(i, wall[i].id, wall[i].point);  
    }
    global_wall_tree.construct();
}

void _set_nd2wall_lvl_mp(std::unordered_map<int, int>& _mp,
                         const std::vector<Node>& all_nds,
                         const GridBTree<int, Point<double>>& global_wall_tree,
                         double D)
//根据节点到物面的距离给节点分等级，其中距离物面D的设为1，距离物面5D的设为2，其他设为3
//输入：_mp：空间节点的id和对应等级的映射
//      all_nds：空间节点坐标
//      global_wall_tree：物面节点的二叉树
//      D：D是五倍最大的位移变形量
{
    int key_temp = 0;       
    double min_d = 0.0;     

    for (const auto& ind : all_nds) 
    {
        //查找ind到物面的最近距离
        global_wall_tree.search(min_d, key_temp, ind.point);
        int lvl = 3;   // 默认等级为 3
        if (min_d <= 5.0 * D) lvl = 2;
        if (min_d <= D)       lvl = 1;

        _mp[ind.id] = lvl;  //把节点 id -> 等级存入 map
    }
}
