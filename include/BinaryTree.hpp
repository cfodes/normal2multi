/* **************************************************************************************************** 
 * Author : –ª¡¡
 * First Created : 2012-10-04,15:53
 * Last Modify   :
 *
 * Binary Tree implemented by template teachnique.
 * ***************************************************************************************************** 
 * */

#pragma once

#include <boost/smart_ptr.hpp>


// class of node in binary tree.
template<typename _Ty>
class BinaryNode
{
	_Ty m_element;			// datas contained.

public:

	BinaryNode<_Ty> *m_left;	// pointer points to left child
	BinaryNode<_Ty> *m_right;	// pointer points to right child

	BinaryNode()
	{
		m_left = NULL;
		m_right = NULL;
	}

    inline bool isLeafNode() const
    {
        return (NULL==m_left && NULL==m_right);
    }

    inline BinaryNode<_Ty> const *left() const
    {
        return m_left;
    }
    inline BinaryNode<_Ty> * left() 
    {
        return m_left;
    }
    inline _Ty const &leftElem() const
    {
        return m_left->m_element;
    }
    inline _Ty &leftElem() 
    {
        return m_left->m_element;
    }
    inline BinaryNode<_Ty> const *right() const
    {
        return m_right;
    }
    inline BinaryNode<_Ty> *right() 
    {
        return m_right;
    }
    inline _Ty const &rightElem() const
    {
        return m_right->m_element;
    }
    inline _Ty &rightElem() 
    {
        return m_right->m_element;
    }
    inline _Ty const &elem() const
    {
        return m_element;
    }
    inline _Ty &elem() 
    {
        return m_element;
    }
};


template<typename _Ty>
class FullBinaryTree
{
public:
    typedef BinaryNode<_Ty> node_type;
    typedef BinaryNode<_Ty> * node_ptr;
	node_type *m_root;	// root node

    inline node_type * root()
    {
        return m_root;
    }
    inline node_type const * root() const
    {
        return m_root;
    }

    inline bool empty() const
    {
        return (NULL==m_root);
    }

public:
	FullBinaryTree()
	{
		m_root = NULL;
	}
	~FullBinaryTree()
	{
		deallocate(m_root);
	}

    void deallocate()
    {
        if (this->empty())
        {
            return;
        }
        deallocate(m_root);
        delete m_root;
        m_root = NULL;
    }

	void deallocate(node_type *node)
	{
		if (NULL == node)
		{
			return ;
		}
		if (NULL != node->m_left)
		{
			deallocate(node->m_left);
			delete node->m_left;
            node->m_left = NULL;
		}
		if (NULL != node->m_right)
		{
			deallocate(node->m_right);
			delete node->m_right;
            node->m_right = NULL;
		}
	}
};

template<typename _Ty>
class BNode
{
    _Ty m_element;

    public:

    boost::shared_ptr<BNode<_Ty> > m_lptr;
    boost::shared_ptr<BNode<_Ty> > m_rptr;

    BNode(_Ty const &item = _Ty(), 
            boost::shared_ptr<BNode<_Ty> > lptr =  boost::shared_ptr<BNode<_Ty> >(),
            boost::shared_ptr<BNode<_Ty> > rptr =  boost::shared_ptr<BNode<_Ty> >()):
        m_element(item), m_lptr(lptr), m_rptr(rptr) {}

    inline bool isLeafNode() const
    {
        return (NULL==m_lptr.get() && NULL==m_rptr.get());
    }

    inline _Ty const &leftElem() const
    {
        return m_lptr->m_element;
    }
    inline _Ty &leftElem() 
    {
        return m_lptr->m_element;
    }
    inline _Ty const &rightElem() const
    {
        return m_rptr->m_element;
    }
    inline _Ty &rightElem() 
    {
        return m_rptr->m_element;
    }
    inline _Ty const &elem() const
    {
        return m_element;
    }
    inline _Ty &elem() 
    {
        return m_element;
    }

    inline boost::shared_ptr<BNode<_Ty> > const &lptr() const
    {
        return m_lptr;
    }
    inline boost::shared_ptr<BNode<_Ty> > &lptr() 
    {
        return m_lptr;
    }
    inline boost::shared_ptr<BNode<_Ty> > const &rptr() const
    {
        return m_rptr;
    }
    inline boost::shared_ptr<BNode<_Ty> > &rptr() 
    {
        return m_rptr;
    }
};

template<typename _Ty>
class BTree
{
    boost::shared_ptr<BNode<_Ty> > m_root;	// root node

    public:

    typedef BNode<_Ty> node_type;
    typedef boost::shared_ptr<BNode<_Ty> > node_ptr;

    BTree() : m_root(boost::shared_ptr<BNode<_Ty> >())
    {
    }

    void deallocate()
    {
        if (this->empty())
        {
            return;
        }
        deallocate(m_root);
        m_root.reset();
    }

	void deallocate(node_ptr &node)
	{
		if (NULL==node.get())
		{
			return ;
		}
		if (NULL!=node->m_lptr.get())
		{
			deallocate(node->m_lptr);
			node->m_lptr.reset();
		}
		if (NULL!=node->m_rptr.get())
		{
			deallocate(node->m_rptr);
			node->m_rptr.reset();
		}
	}

    inline bool empty() const
    {
        return (NULL==m_root.get());
    }

    inline boost::shared_ptr<BNode<_Ty> > &root()
    {
        return m_root;
    }
    inline boost::shared_ptr<BNode<_Ty> > const &root() const
    {
        return m_root;
    }
};


