/*
 Copyright 2013 Edouard Griffiths <f4exb at free dot fr>

 This file is part of CCSoft. A Convolutional Codes Soft Decoding library

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Boston, MA  02110-1301  USA

 Node in the convolutional code tree

 */

#ifndef __CC_TREE_NODE_H__
#define __CC_TREE_NODE_H__

#include <vector>

namespace ccsoft
{

template<typename T_IOSymbol, typename T_Register, typename T_EdgeTag>
class CC_TreeEdge;

/**
 * \brief A node in the code tree
 * \tparam T_IOSymbol Type of the input and output symbols
 * \tparam T_Register Type of the encoder internal registers
 */
template<typename T_IOSymbol, typename T_Register, typename T_EdgeTag>
class CC_TreeNode
{
public:
    /**
     * Constructor
     * \param _id Unique ID of the edge
     * \param _p_incoming_edge Pointer to the incoming edge to the node
     * \param _path_metric Path metric at the node
     * \param _depth This node depth
     */
    CC_TreeNode(unsigned int _id,
            CC_TreeEdge<T_IOSymbol, T_Register, T_EdgeTag> *_p_incoming_edge,
            float _path_metric,
            int _depth) :
                id(_id),
                p_incoming_edge(_p_incoming_edge),
                path_metric(_path_metric),
                depth(_depth),
                on_final_path(false)
    {}

    /**
     * Destructor. Destroys all outgoing edges
     */
    ~CC_TreeNode()
    {
        delete_outgoing_edges();
    }

    /**
     * Add an outgoing edge
     */
    void add_outgoing_edge(CC_TreeEdge<T_IOSymbol, T_Register, T_EdgeTag> *p_outgoing_edge)
    {
        p_outgoing_edges.push_back(p_outgoing_edge);
    }

    /**
     * Delete outgoing edges
     */
    void delete_outgoing_edges()
    {
        typename std::vector<CC_TreeEdge<T_IOSymbol, T_Register, T_EdgeTag>*>::iterator e_it = p_outgoing_edges.begin();
    
        for (; e_it != p_outgoing_edges.end(); ++e_it)
        {
            if (*e_it)
            {
                delete *e_it;
                *e_it = 0;
            }
        }

        p_outgoing_edges.clear();
    }

    /**
     * Return a R/O reference to the outgoing edges
     */
    const std::vector<CC_TreeEdge<T_IOSymbol, T_Register, T_EdgeTag>*>& get_outgoing_edges() const
    {
        return p_outgoing_edges;
    }
    
    /**
     * Return a R/W reference to the outgoing edges
     */
    std::vector<CC_TreeEdge<T_IOSymbol, T_Register, T_EdgeTag>*>& get_outgoing_edges() 
    {
        return p_outgoing_edges;
    }
    
    /**
     * Get pointer to the incoming edge
     */
    CC_TreeEdge<T_IOSymbol, T_Register, T_EdgeTag> *get_incoming_edge()
    {
        return p_incoming_edge;
    }

    /**
     * Get path metric to the node
     */
    float get_path_metric() const
    {
        return path_metric;
    }

    /**
     * Get the depth of the node
     */
    int get_depth() const
    {
        return depth;
    }
    
    /**
     * Get node id
     */
    unsigned int get_id() const
    {
        return id;
    }

    /**
     * For ordering by increasing path metric
     */
    bool operator<(const CC_TreeNode<T_IOSymbol, T_Register, T_EdgeTag>& other) const
    {
        return path_metric < other.path_metric;
    }

    /**
     * For ordering by decreasing path metric
     */
    bool operator>(const CC_TreeNode<T_IOSymbol, T_Register, T_EdgeTag>& other) const
    {
        return path_metric > other.path_metric;
    }

    /**
     * Get saved encoder registers reference
     */
    const std::vector<T_Register>& get_registers() const
    {
        return registers;
    }

    /**
     * Save encoder registers
     */
    void set_registers(const std::vector<T_Register>& _registers)
    {
        registers = _registers;
    }
    
    /**
     * Set the "on final path" marker
     */
    void set_on_final_path(bool _on_final_path = true)
    {
        on_final_path = _on_final_path;
    }
    
    /**
     * Test the "on final path" marker
     */
    bool is_on_final_path()
    {
        return on_final_path;
    }
    
    /**
     * Ordering - lesser
     */
    bool operator<(const CC_TreeNode<T_IOSymbol, T_Register, T_EdgeTag>& other)
    {
        if (path_metric == other.path_metric)
        {
            return id < other.id;
        }
        else
        {
            return path_metric < other.path_metric;
        }
    }

    /**
     * Ordering - greater
     */
    bool operator>(const CC_TreeNode<T_IOSymbol, T_Register, T_EdgeTag>& other)
    {
        if (path_metric == other.path_metric)
        {
            return id > other.id;
        }
        else
        {
            return path_metric > other.path_metric;
        }
    }

protected:
    unsigned int id; //!< Node's unique ID
    std::vector<CC_TreeEdge<T_IOSymbol, T_Register, T_EdgeTag>*> p_outgoing_edges; //!< Outgoing edges pointers
    CC_TreeEdge<T_IOSymbol, T_Register, T_EdgeTag> *p_incoming_edge; //!< Pointer to the incoming edge
    float path_metric; //!< Path metric to the node
    int depth; //!< Depth of node in the tree: 0 = root
    std::vector<T_Register> registers; //!< state of encoder registers at node
    bool on_final_path; //!< Marks node when backtracking the solution
};

} // namespace ccsoft

#endif // __CC_TREE_NODE_H__
