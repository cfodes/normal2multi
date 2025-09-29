#pragma once

//#include <Gnl/src/vector.h>
//#include <tools.hpp>
//#include <gnl/src/point.h>

#include <limits>
#include <utility>

#include "BinaryTree.hpp"

template<typename _PTTy>
class Block
{
    int m_id0;
    int m_idM;

    _PTTy m_max;
    _PTTy m_min;

    public:

    Block() :
        m_id0(0),
        m_idM(0)
    {
    }

    inline void setIndex(int const &id0, int const &idM)
    {
        m_id0 = id0;
        m_idM = idM;
    }

    inline int const &id0() const
    {
        return m_id0;
    }
    inline int const &idM() const
    {
        return m_idM;
    }

    inline _PTTy const &max() const
    {
        return m_max;
    }
    inline _PTTy &max() 
    {
        return m_max;
    }
    inline _PTTy const &min() const
    {
        return m_min;
    }
    inline _PTTy &min() 
    {
        return m_min;
    }
};


namespace gnl
{
    // 求解两个变量中的最大值
    template<typename _T1, typename _T2>
    inline auto Max(const _T1& a, const _T2& b)
    {
        return a > b ? a : b;
    }

    // 求解两个变量中的最小值
    template<typename _T1, typename _T2>
    inline auto Min(const _T1& a, const _T2& b)
    {
        return a > b ? b : a;
    }

    // 求解绝对值
    template<typename _Ty>
    inline _Ty ABS(const _Ty& a)
    {
        return a > _Ty() ? a : -a;
    }
}

// ===================================================================================================================
// Binary Tree for the searching of the minimal distance between a point and all of the point on the wall.
// ___________________________________________________________________________________________________________________
// template parameters:
//  _KEY  : a key to identity the points on the wall.
//          . overload the operator =.
//  _PTTy : type of point. should provide:
//          . provide ::value_type, which is the type of each component of the point.
//          . overload the copy constructor.
//          . overload the operator =, -=.
//          . overload the operator[] to access each components of the points.
//          . provide size() function to access the dimension of the points.
//          . provide distSquare(_PTTy p1, _PTTy p2) to compute the square of the distance between p1 and p2.
//          . provide square(_PTTy p1) to compute the sume of square of each component of p1.
// ___________________________________________________________________________________________________________________
template<typename _KEY, typename _PTTy>
class GridBTree : private BTree<Block<_PTTy> >
{
    typedef typename _PTTy::value_type _REAL_T;
    typedef typename BTree<Block<_PTTy> >::node_type node_type;
    typedef typename BTree<Block<_PTTy> >::node_ptr node_ptr;

#ifdef _DEBUG
    private : int m_count;
#endif

    public: GridBTree()
            {
#ifdef _DEBUG
                m_count = 0;
#endif
            }

    public : void setGridSize(int const &size)
             {
                 m_grids.resize(size);
             }
    public : void setGrid(int const &id, _KEY const &key, _PTTy const &grid)
             {
                 m_grids[id].first = key;
                 m_grids[id].second = grid;
             }

    public : void construct()
             {
                 if (!this->empty())
                 {
                     printf("EMPTY\n");
                     this->deallocate();
                 }
                 //printf("START SORT\n");
                 sort(this->root());
                 //printf("SET BOX\n");
                 setBox(*(this->root()));
             }

             // ========================================================================================
             // search minimal distance between the point and the grids on the wall.
             // ________________________________________________________________________________________
             // pars:
             // 	d 	  :[o] minimal distance.
             // 	key   :[o] the index of the point on wall, the distance between which and point is minimal.
             // 	point :[i] point we want to compute the distance between which and wall.
             // return
             //  VOID
             // ________________________________________________________________________________________
    public : void search(_REAL_T &d, _KEY &key, _PTTy const &point) const
             {
                 _REAL_T dd = (std::numeric_limits<_REAL_T>::max)();
                 this->search(dd, key, point, *(this->root()));
                 d = sqrt(dd);
             }

    public : inline int grid_size() const
             {
                 return m_grids.size();
             }

    public : inline _PTTy const &grid(int const &id) const
             {
                 return m_grids[id].second;
             }

    public : inline _KEY const &key(int const &id) const
             {
                 return m_grids[id].first;
             }

    private:

             std::vector<std::pair<_KEY, _PTTy>> m_grids;

             // ========================================================================================
             // search minimal distance between the point and the grids on the wall.
             // ________________________________________________________________________________________
             // pars:
             // 	dd	  :[o] square of minimal distance.
             // 	key   :[o] the index of the point on wall, the distance between which and point is minimal.
             // 	point :[i] point we want to compute the distance between which and wall.
             // 	node  :[i] node we start the searching procedure.
             // return 
             //  VOID
             // ________________________________________________________________________________________
    private : void search(_REAL_T &dd, _KEY &key, _PTTy const &point, node_type const &node) const
              {
                  if (node.isLeafNode())
                  {
                      _REAL_T ddTmp;
                      for (int id=node.elem().id0(); id<node.elem().idM(); ++id)
                      {
                          //ddTmp = distSquare(m_grids[id].second, point);
                          _REAL_T Re = _REAL_T();
                          Re += pow(m_grids[id].second.x - point.x, 2);
                          Re += pow(m_grids[id].second.y - point.y, 2);
                          Re += pow(m_grids[id].second.z - point.z, 2);
                          ddTmp = Re;

                          if (ddTmp<dd)
                          {
                              dd = ddTmp;
                              key = m_grids[id].first;
                          }
                      }
                  }
                  else
                  {
                      _PTTy delta;

                      // compute distance to left block.
                      _PTTy dR(point); 
                      _PTTy dL(point); 
                      dR -= node.leftElem().max();
                      dL -= node.leftElem().min();
                      for (int i=0; i!=_PTTy().size(); ++i)
                      {
                          delta[i] = (dL[i] * dR[i] < 0) ? 0 : gnl::Min(gnl::ABS(dL[i]), gnl::ABS(dR[i]));
                      }
                      _REAL_T const ddL = square(delta);

                      // compute distance to right block.
                      dR = point; 
                      dL = point; 
                      dR -= node.rightElem().max();
                      dL -= node.rightElem().min();
                      for (int i=0; i!=_PTTy().size(); ++i)
                      {
                          delta[i] = (dL[i] * dR[i] < 0) ? 0 : gnl::Min(gnl::ABS(dL[i]), gnl::ABS(dR[i]));
                      }
                      _REAL_T const ddR = square(delta);

                      // compare.
                      if (dd<ddL && dd<ddR)
                      {
                          return ;
                      }
                      if (ddL<ddR)
                      {
                          search(dd, key, point, *node.lptr());
                          if (dd<ddR)
                          {
                              return ;
                          }
                          search(dd, key, point, *node.rptr());
                      }
                      else
                      {
                          search(dd, key, point, *node.rptr());
                          if (dd<ddL)
                          {
                              return ;
                          }
                          search(dd, key, point, *node.lptr());
                      }
                  }
              }

              // ========================================================================================
              // sort grids by the fast sort method.
              // ________________________________________________________________________________________
              // pars:
              // 	pNode	:[i] The node of which we start sorting.
              // return
              //  VOID
              // ________________________________________________________________________________________
    private : void sort(node_ptr &pNode)
              {
#ifdef _DEBUG
                  ++m_count;
                  //printf("LINE:%d, %d\n", __LINE__, m_count);
#endif
                  if (this->empty())
                  {
                      // printf("LINE:%d\n", __LINE__);
                      this->root() = boost::make_shared<node_type>();
                      // printf("LINE:%d\n", __LINE__);
                      this->root()->elem().setIndex(0, m_grids.size());
                      // printf("LINE:%d\n", __LINE__);
                      sort(this->root());
                      // printf("LINE:%d\n", __LINE__);
                      return ;
                  }

                  int const &id0 = pNode->elem().id0();
                  int const &idM = pNode->elem().idM();
                  // error ! when the value is smaller than 3, an debug will occur.
                  if (idM-id0<=5)
                  {
                      return ;
                  }
                  // set box for all blocks.
                  _PTTy min;
                  _PTTy max;
                  _PTTy mid;
                  _PTTy const *pPnt;

                  // printf("LINE:%d\n", __LINE__);
                  min = m_grids[id0].second;
                  // printf("LINE:%d\n", __LINE__);
                  max = min;
                  // printf("LINE:%d\n", __LINE__);
                  for (int id=id0; id<idM; ++id)
                  {
                      // printf("LINE:%d\n", __LINE__);
                      pPnt = &(m_grids[id].second);

                      // printf("LINE:%d\n", __LINE__);
                      for (int i=0; i!=_PTTy().size(); ++i)
                      {
                          // printf("LINE:%d\n", __LINE__);
                          min[i] = gnl::Min(min[i], (*pPnt)[i]);
                          // printf("LINE:%d\n", __LINE__);
                          max[i] = gnl::Max(max[i], (*pPnt)[i]);
                          // printf("LINE:%d\n", __LINE__);
                      }
                  }
                  int dim = 0;
                  _REAL_T dimMax = 0;
                  // printf("LINE:%d\n", __LINE__);
                  for (int i=0; i<_PTTy().size(); ++i)
                  {
                      // printf("LINE:%d\n", __LINE__);
                      _REAL_T const space = max[i] - min[i];
                      // printf("LINE:%d\n", __LINE__);
                      if (space>dimMax)
                      {
                          // printf("LINE:%d\n", __LINE__);
                          dimMax = space;
                          // printf("LINE:%d\n", __LINE__);
                          dim = i;
                      }
                  }
                  // printf("LINE:%d\n", __LINE__);
                  mid[dim] = (max[dim] + min[dim]) / 2;

                  // sorting by the fast sort method.
                  int i = id0;
                  int j = idM - 1;
                  do 
                  {
                      // printf("LINE:%d\n", __LINE__);
                      for (; i<idM; ++i)
                      {
                          // printf("LINE:%d\n", __LINE__);
                          if (m_grids[i].second[dim] > mid[dim])
                          {
                              // printf("LINE:%d\n", __LINE__);
                              break;
                          }
                      }
                      // printf("LINE:%d\n", __LINE__);
                      for (; j>=id0; --j)
                      {
                          // printf("LINE:%d\n", __LINE__);
                          if (m_grids[j].second[dim] < mid[dim])
                          {
                              // printf("LINE:%d\n", __LINE__);
                              break;
                          }
                      }
                      // printf("LINE:%d\n", __LINE__);
                      if (i<j)
                      {
                          // printf("LINE:%d\n", __LINE__);
                          std::swap(m_grids[i], m_grids[j]);
                      }
                      else
                      {
                          break;
                      }
                  }
                  while (true);
                  // printf("LINE:%d\n", __LINE__);
                  std::swap(m_grids[i], m_grids[idM-1]);

                  // construct left child node.
                  // printf("LINE:%d\n", __LINE__);
                  pNode->lptr() = boost::make_shared<node_type>();
                  // printf("LINE:%d\n", __LINE__);
                  pNode->leftElem().setIndex(id0, i);
                  // printf("LINE:%d\n", __LINE__);
                  // printf("LPTR:%d\n", pNode->lptr().get());
                  sort(pNode->lptr());

                  // construct right child node.
                  // printf("LINE:%d\n", __LINE__);
                  pNode->rptr() = boost::make_shared<node_type>();
                  // printf("LINE:%d\n", __LINE__);
                  pNode->rightElem().setIndex(i, idM);
                  // printf("LINE:%d\n", __LINE__);
                  sort(pNode->rptr());
                  // printf("LINE:%d\n", __LINE__);
              }

              // ========================================================================================
              //  compute the box cover each blocks associated with the node.
              // ________________________________________________________________________________________
              //  pars:
              //  	node :[i] the node is associated with the block the box of which will be computed.
              // ________________________________________________________________________________________
    private : void setBox(node_type &node)
              {
                  if (node.isLeafNode())
                  {
                      int const &id0 = node.elem().id0();
                      int const &idM = node.elem().idM();
                      _PTTy &max = node.elem().max();
                      _PTTy &min = node.elem().min();
                      min = max = this->m_grids[id0].second;
                      for (int id=id0; id<idM; ++id)
                      {
                          _PTTy const &grd = m_grids[id].second;
                          for (int i=0; i!=_PTTy().size(); ++i)
                          {
                              min[i] = gnl::Min(min[i], grd[i]);
                              max[i] = gnl::Max(max[i], grd[i]);
                          }
                      }
                  }
                  else
                  {
                      setBox(*node.lptr());
                      setBox(*node.rptr());
                      _PTTy &max = node.elem().max();
                      _PTTy &min = node.elem().min();
                      for (int i=0; i!=_PTTy().size(); ++i)
                      {
                          max[i] = gnl::Max(node.leftElem().max()[i], node.rightElem().max()[i]);
                          min[i] = gnl::Min(node.leftElem().min()[i], node.rightElem().min()[i]);
                      }
                  }
              }
};
