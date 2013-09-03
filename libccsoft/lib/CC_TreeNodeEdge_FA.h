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

 Combination of a Node and its incoming Edge in the convolutional code tree
 (node+edge combo). In a tree structure you don't need to store nodes and
 edges as nodes have a single incoming edge. So a node can incorporate
 its incoming edge.

 This version uses fixed arrays for registers and forward node-edges

 */
#ifndef __CC_TREE_NODE_EDGE_FA_H__
#define __CC_TREE_NODE_EDGE_FA_H__

#include "CC_TreeNodeEdge_base.h"
#include "CC_EncodingRegisters_FA.h"

#include <vector>
#include <array>
#include <algorithm>

namespace ccsoft
{

/**
 * \brief Represents a node and its incoming edge in the code tree
 * This version uses a fixed array to store forward node+edges pointers.
 * N_k template parameter gives the size of the input symbol (k parameter).
 * There are (1<<N_k) forward node+edges.
 * \tparam T_IOSymbol Type of the input and output symbols
 * \tparam T_Register Type of the encoder internal registers
 * \tparam T_Tag Type of the node-edge tag
 * \tparam N_k Input symbol size in bits (k parameter)
 */
template<typename T_IOSymbol, typename T_Register, typename T_Tag, unsigned int N_k>
class CC_TreeNodeEdge_FA : public CC_TreeNodeEdge_base<T_IOSymbol, T_Tag>, public CC_EncodingRegisters_FA<T_Register, N_k>
{

public:
    /**
     * Constructor
     * \param _id Unique ID of the edge
     * \param _p_incoming_edge Pointer to the incoming edge to the node
     * \param _in_symbol Input symbol corresponding to the edge
     * \param _metric Metric of the edge
     * \param _path_metric Path metric at the node
     * \param _depth This node depth
     */
	CC_TreeNodeEdge_FA(unsigned int _id,
			CC_TreeNodeEdge_FA<T_IOSymbol, T_Register, T_Tag, N_k> *_p_incoming_node_edge,
			const T_IOSymbol& _in_symbol,
            float _incoming_edge_metric,
            float _path_metric,
            int _depth) :
                CC_TreeNodeEdge_base<T_IOSymbol, T_Tag>(_id, _in_symbol, _incoming_edge_metric, _path_metric, _depth),
                p_incoming_node_edge(_p_incoming_node_edge)
    {
        clear_outgoing_edges();
    }

	/**
	 * Destructor
	 */
	~CC_TreeNodeEdge_FA()
	{
		// deletes all outgoing edge+nodes
		delete_outgoing_node_edges();
	}

    /**
     * Add an outgoing edge
     * \param p_outgoing_node_edge Outgoing edge+node
     */
    void set_outgoing_node_edge(CC_TreeNodeEdge_FA<T_IOSymbol, T_Register, T_Tag, N_k> *p_outgoing_node_edge, unsigned int index)
    {
    	p_outgoing_node_edges[index] = p_outgoing_node_edge;
    }

    /**
     * Delete outgoing edges
     */
    void delete_outgoing_node_edges()
    {
        typename std::array<CC_TreeNodeEdge_FA<T_IOSymbol, T_Register, T_Tag, N_k>*, (1<<N_k)>::iterator ne_it = p_outgoing_node_edges.begin();

        for (; ne_it != p_outgoing_node_edges.end(); ++ne_it)
        {
            if (*ne_it)
            {
                delete *ne_it;
                *ne_it = 0;
            }
        }

        clear_outgoing_edges();
    }

    /**
     * Verifies validity of outgoing node+edges i.e. pointers are all non null
     */
    bool valid_outgoing_node_edges(unsigned int index_limit = (1<<N_k))
    {
        for (unsigned int i=0; i<(1<<N_k); i++)
        {
            if (i == index_limit)
            {
                break;
            }

            if (p_outgoing_node_edges[i] == 0)
            {
                return false;
            }
        }

        return true;
    }

    /**
     * Return a R/O reference to the outgoing node+edges
     */
    const std::array<CC_TreeNodeEdge_FA<T_IOSymbol, T_Register, T_Tag, N_k>*, (1<<N_k)>& get_outgoing_node_edges() const
    {
        return p_outgoing_node_edges;
    }

    /**
     * Return a R/W reference to the outgoing edges
     */
    std::array<CC_TreeNodeEdge_FA<T_IOSymbol, T_Register, T_Tag, N_k>*, (1<<N_k)>& get_outgoing_node_edges()
    {
        return p_outgoing_node_edges;
    }

    /**
     * Get pointer to the incoming edge
     */
    CC_TreeNodeEdge_FA<T_IOSymbol, T_Register, T_Tag, N_k> *get_incoming_node_edge()
    {
        return p_incoming_node_edge;
    }
    
protected:
    void clear_outgoing_edges()
    {
        std::fill(p_outgoing_node_edges.begin(), p_outgoing_node_edges.end(), (CC_TreeNodeEdge_FA<T_IOSymbol, T_Register, T_Tag, N_k>*) 0);
        //p_outgoing_node_edges.fill(0);
    }

    std::array<CC_TreeNodeEdge_FA<T_IOSymbol, T_Register, T_Tag, N_k>*, (1<<N_k)> p_outgoing_node_edges; //!< Outgoing edges+node pointers
    CC_TreeNodeEdge_FA<T_IOSymbol, T_Register, T_Tag, N_k> *p_incoming_node_edge; //!< Pointer to the incoming edge+node
};

} // namespace ccsoft

#endif
