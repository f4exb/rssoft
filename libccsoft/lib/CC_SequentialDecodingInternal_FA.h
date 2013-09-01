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

 Convolutional soft-decision sequential decoder generic (virtual) class
 Based on node+edge combination in the code tree

 */
#ifndef __CC_SEQUENTIAL_DECODING_INERNAL_FA_H__
#define __CC_SEQUENTIAL_DECODING_INERNAL_FA_H__

#include "CC_TreeNodeEdge_FA.h"
#include "CC_ReliabilityMatrix.h"
#include "CC_TreeGraphviz_FA.h"

#include <cmath>
#include <algorithm>
#include <iostream>



namespace ccsoft
{

/**
 * \brief Convolutional soft-decision sequential decoder generic (virtual) class for algorithm internal use.
 * It is tainted by the type of code tree node+edge tag that is algorithm dependant. It contains the code tree root node and some
 * common methods.
 * This version uses a fixed array to store forward node+edges pointers.
 * N_k template parameter gives the size of the input symbol (k parameter).
 * There are (1<<N_k) forward node+edges.
 * \tparam T_Register Type of the encoder internal registers
 * \tparam T_IOSymbol Type of the input and output symbols
 * \tparam T_EdgeTag Type of the code tree node+edge tag
 * \tparam N_k Size of an input symbol in bits (k parameter)
 */
template<typename T_Register, typename T_IOSymbol, typename T_Tag, unsigned int N_k>
class CC_SequentialDecodingInternal_FA
{
public:
    /**
     * Constructor
     * \param constraints Vector of register lengths (constraint length + 1). The number of elements determines k.
     * \param genpoly_representations Generator polynomial numeric representations. There are as many elements as there
     * are input bits (k). Each element is itself a vector with one polynomial value per output bit. The smallest size of
     * these vectors is retained as the number of output bits n. The input bits of a symbol are clocked simultaneously into
     * the right hand side, or least significant position of the internal registers. Therefore the given polynomial representation
     * of generators should follow the same convention.
     */
	CC_SequentialDecodingInternal_FA() :
        root_node(0)
	{}

	/**
	 * Destructor
	 */
	virtual ~CC_SequentialDecodingInternal_FA()
	{
        if (root_node)
        {
            delete root_node;
            root_node = 0;
        }
	}

    /**
     * Reset the decoding process
     */
    void reset()
    {
        if (root_node)
        {
            delete root_node;
            root_node = 0;
        }
    }

protected:
    /**
     * Base 2 logarithm
     */
    float log2(float x)
    {
        return log(x)/log(2.0);
    }

    /**
     * Initialize process at the root node
     */
    void init_root()
    {
        root_node = new CC_TreeNodeEdge_FA<T_IOSymbol, T_Register, T_Tag, N_k>(0, 0, 0, 0.0, 0.0, -1);
    }

    /**
     * Visit a new node in the code tree
     * \param node Node to visit
     * \param relmat Reliability matrix reference
     */
    virtual void visit_node_forward(CC_TreeNodeEdge_FA<T_IOSymbol, T_Register, T_Tag, N_k>* node, const CC_ReliabilityMatrix& relmat) = 0;

    /**
     * Back track from a node. When the node is the selected terminal node it is used to retrieve the decoded message
     * \param node Node to track back from
     * \param decoded_message Symbols corresponding to the edge ordered from root node to the given node
     * \param mark_nodes Mark the nodes along the path
     */
    void back_track(CC_TreeNodeEdge_FA<T_IOSymbol, T_Register, T_Tag, N_k>* node_edge, std::vector<T_IOSymbol>& decoded_message, bool mark_nodes = false)
    {
        std::vector<T_IOSymbol> reversed_message;
        CC_TreeNodeEdge_FA<T_IOSymbol, T_Register, T_Tag, N_k> *cur_node_edge = node_edge;
        CC_TreeNodeEdge_FA<T_IOSymbol, T_Register, T_Tag, N_k> *incoming_node_edge;

        reversed_message.push_back(cur_node_edge->get_in_symbol());

        while (incoming_node_edge = (cur_node_edge->get_incoming_node_edge()))
        {
            cur_node_edge->set_on_final_path(mark_nodes);

            if (incoming_node_edge->get_depth() >= 0) // don't take root node
            {
                reversed_message.push_back(incoming_node_edge->get_in_symbol());
            }

            cur_node_edge = incoming_node_edge;
        }

        decoded_message.resize(reversed_message.size());
        std::reverse_copy(reversed_message.begin(), reversed_message.end(), decoded_message.begin());
    }

    /**
     * Print the code tree in Graphviz dot format
     * \param os Output stream
     */
    void print_dot_internal(std::ostream& os)
    {
        CC_TreeGraphviz_FA<T_IOSymbol, T_Register, T_Tag, N_k>::create_dot(root_node, os);
    }
    
    CC_TreeNodeEdge_FA<T_IOSymbol, T_Register, T_Tag, N_k> *root_node; //!< Root node
};


} // namespace ccsoft


#endif // __CC_SEQUENTIAL_DECODING_INERNAL_H__
