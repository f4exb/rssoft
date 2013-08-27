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

 */
#ifndef __CC_SEQUENTIAL_DECODING_H__
#define __CC_SEQUENTIAL_DECODING_H__

#include "CC_Encoding.h"
#include "CC_TreeEdge.h"
#include "CC_TreeNode.h"
#include "CC_ReliabilityMatrix.h"
#include "CC_TreeGraphviz.h"

#include <cmath>
#include <algorithm>
#include <iostream>



namespace ccsoft
{

/**
 * Base 2 logarithm
 */
float log2(float x)
{
    return log(x)/log(2.0);
}

/**
 * \brief class used for node ordering
 */
class NodeOrdering
{
public:
    NodeOrdering(float _path_metric, unsigned int _node_id) :
        path_metric(_path_metric),
        node_id(_node_id)
    {}

    ~NodeOrdering()
    {}

    bool operator>(const NodeOrdering& other) const
    {
        if (path_metric == other.path_metric)
        {
            return node_id > other.node_id;
        }
        else
        {
            return path_metric > other.path_metric;
        }
    }

    float path_metric;
    unsigned int node_id;
};

template<typename T_Node>
bool node_pointer_ordering(T_Node* n1, T_Node* n2)
{
    if (n1->get_path_metric() == n2->get_path_metric())
    {
        return n1->get_id() > n2->get_id();
    }
    else
    {
        return n1->get_path_metric() > n2->get_path_metric();
    }
}

/**
 * \brief Convolutional soft-decision sequential decoder generic (virtual) class. This is the public interface.
 * \tparam T_Register Type of the encoder internal registers
 * \tparam T_IOSymbol Type of the input and output symbols
 */
template<typename T_Register, typename T_IOSymbol>
class CC_SequentialDecoding
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
	CC_SequentialDecoding(const std::vector<unsigned int>& constraints,
            const std::vector<std::vector<T_Register> >& genpoly_representations) :
                encoding(constraints, genpoly_representations),
                use_metric_limit(false),
                metric_limit(0.0),
                use_node_limit(false),
                node_limit(0),
                codeword_score(0.0),
                cur_depth(-1),
                max_depth(0),
                node_count(0),
                edge_count(0),
                tail_zeros(true),
                edge_bias(0.0),
                verbosity(0)
	{}

	/**
	 * Destructor
	 */
	virtual ~CC_SequentialDecoding()
	{
	}

    /**
     * Set the node limit threshold
     */
    void set_node_limit(unsigned int _node_limit)
    {
        node_limit = _node_limit;
        use_node_limit = true;
    }

    /**
     * Reset the node limit threshold. The process will continue until out of memory or end of the tree.
     */
    void reset_node_limit()
    {
        use_node_limit = false;
    }

    /**
     * Set the metric limit threshold
     */
    void set_metric_limit(float _metric_limit)
    {
        metric_limit = _metric_limit;
        use_metric_limit = true;
    }

    /**
     * Reset the metric limit threshold. The process will continue until out of memory or end of the tree.
     */
    void reset_metric_limit()
    {
        use_metric_limit = false;
    }

    /**
     * Set the tail zeros option
     */
    void set_tail_zeros(bool _tail_zeros)
    {
        tail_zeros = _tail_zeros;
    }

    /**
     * Reset the decoding process
     */
    void reset()
    {
        node_count = 0;
        edge_count = 0;
        codeword_score = 0.0;
        cur_depth = -1;
        max_depth = 0;
        encoding.clear(); // clear encoder's registers
    }

    /**
     * Get encoding object reference
     */
    CC_Encoding<T_Register, T_IOSymbol>& get_encoding()
    {
        return encoding;
    }

    /**
     * Get the codeword score. Valid only if decode returned successfully.
     */
    float get_score() const
    {
        return codeword_score;
    }

    /**
     * Get the codeword score in dB/Symbol units. Valid only if decode returned successfully.
     */
    float get_score_db_sym() const
    {
        if (cur_depth > 0)
        {
            return (10.0*log(2.0)*codeword_score) / cur_depth;
        }
        else
        {
            return 0.0;
        }
    }

    /**
     * Get the number of nodes created minus root node
     */
    unsigned int get_nb_nodes() const
    {
        return node_count;
    }

    /**
     * Get the current depth
     */
    unsigned int get_current_depth() const
    {
    	return cur_depth;
    }

    /**
     * Get the maximum depth reached
     */
    unsigned int get_max_depth() const
    {
    	return max_depth;
    }
    
    /**
     * Set the edge metric bias
     */
    void set_edge_bias(float _edge_bias)
    {
        edge_bias = _edge_bias;
    }

    /**
     * Set verbosity level
     */
    void set_verbosity(unsigned int _verbosity)
    {
        verbosity = _verbosity;
    }

    /**
     * Print the dot (Graphviz) file of the current decode tree to an output stream
     * \param os Output stream
     */
    virtual void print_dot(std::ostream& os) = 0;
    
    /**
     * Print statistics to an output stream
     * \param os Output stream
     * \param success true if decoding was successful
     */
    virtual void print_stats(std::ostream& os, bool success) = 0;

    /**
     * Decodes given the reliability matrix
     * \param relmat Reference to the reliability matrix
     * \param decoded_message Vector of symbols of retrieved message
     */
    virtual bool decode(const CC_ReliabilityMatrix& relmat, std::vector<T_IOSymbol>& decoded_message) = 0;

protected:
    CC_Encoding<T_Register, T_IOSymbol> encoding;   //!< Convolutional encoding object
    bool use_metric_limit;    //!< True if a give up path metric threshold is used
    float metric_limit;       //!< The give up path metric threshold
    bool use_node_limit;      //!< Stop above number of nodes threshold
    unsigned int node_limit;  //!< Number of nodes threshold
    float codeword_score;     //!< Metric of the codeword found if any
    int cur_depth;            //!< Current depth for the encoder
    int max_depth;            //!< Maximum depth reached in the graph
    unsigned int node_count;  //!< Count of nodes in the code tree
    unsigned int edge_count;  //!< Count of edges in the code tree
    bool tail_zeros;          //!< True if tail of m-1 zeros in the message are assumed. This is the default option.
    float edge_bias;          //!< Edge metric bias subtracted from log2 of reliability of the edge
    unsigned int verbosity;   //!< Verbosity level
};

/**
 * \brief Convolutional soft-decision sequential decoder generic (virtual) class for algorithm internal use.
 * It is tainted by the type of code tree edge tag that is algorithm dependant. It contains the code tree root node and some 
 * common methods.
 * \tparam T_Register Type of the encoder internal registers
 * \tparam T_IOSymbol Type of the input and output symbols
 * \tparam T_EdgeTag Type of the code tree edge tag
 */
template<typename T_Register, typename T_IOSymbol, typename T_EdgeTag>
class CC_SequentialDecodingInternal
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
	CC_SequentialDecodingInternal() :
        root_node(0)
	{}

	/**
	 * Destructor
	 */
	virtual ~CC_SequentialDecodingInternal()
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
     * Initialize process at the root node
     */
    void init_root()
    {
        root_node = new CC_TreeNode<T_IOSymbol, T_Register, T_EdgeTag>(0, 0, 0.0, -1);
    }

    /**
     * Visit a new node in the code tree
     * \param node Node to visit
     * \param relmat Reliability matrix reference
     */
    virtual void visit_node_forward(CC_TreeNode<T_IOSymbol, T_Register, T_EdgeTag>* node, const CC_ReliabilityMatrix& relmat) = 0;

    /**
     * Back track from a node. When the node is the selected terminal node it is used to retrieve the decoded message
     * \param node Node to track back from
     * \param decoded_message Symbols corresponding to the edge ordered from root node to the given node
     * \param mark_nodes Mark the nodes along the path
     */
    void back_track(CC_TreeNode<T_IOSymbol, T_Register, T_EdgeTag>* node, std::vector<T_IOSymbol>& decoded_message, bool mark_nodes = false)
    {
        std::vector<T_IOSymbol> reversed_message;
        CC_TreeNode<T_IOSymbol, T_Register, T_EdgeTag> *cur_node = node;
        CC_TreeEdge<T_IOSymbol, T_Register, T_EdgeTag> *incoming_edge;

        while (incoming_edge = (cur_node->get_incoming_edge()))
        {
            cur_node->set_on_final_path(mark_nodes);
            reversed_message.push_back(incoming_edge->get_in_symbol());
            cur_node = incoming_edge->get_p_origin();
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
        CC_TreeGraphviz<T_IOSymbol, T_Register, T_EdgeTag>::create_dot(root_node, os);
    }
    
    CC_TreeNode<T_IOSymbol, T_Register, T_EdgeTag> *root_node; //!< Root node
};


} // namespace ccsoft


#endif // __CC_SEQUENTIAL_DECODING_H__
