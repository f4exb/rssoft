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
#ifndef __CC_SEQUENTIAL_DECODING_H__
#define __CC_SEQUENTIAL_DECODING_H__

#include "CC_Encoding.h"
#include "CC_ReliabilityMatrix.h"
#include "CC_Interleaver.h"

#include <cmath>
#include <algorithm>
#include <iostream>



namespace ccsoft
{

/**
 * \brief class used for node ordering
 */
class NodeEdgeOrdering
{
public:
	NodeEdgeOrdering(float _path_metric, unsigned int _node_id) :
        path_metric(_path_metric),
        node_id(_node_id)
    {}

    ~NodeEdgeOrdering()
    {}

    bool operator>(const NodeEdgeOrdering& other) const
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

template<typename T_NodeEdge>
bool node_edge_pointer_ordering(T_NodeEdge* n1, T_NodeEdge* n2)
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
class CC_SequentialDecoding : public CC_Interleaver<T_IOSymbol>
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
     * Print statistics summary to an output stream
     * \param os Output stream
     * \param success true if decoding was successful
     */
    virtual void print_stats_summary(std::ostream& os, bool success) = 0;

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
    bool tail_zeros;          //!< True if tail of m-1 zeros in the message are assumed. This is the default option.
    float edge_bias;          //!< Edge metric bias subtracted from log2 of reliability of the edge
    unsigned int verbosity;   //!< Verbosity level
};

} // namespace ccsoft


#endif // __CC_SEQUENTIAL_DECODING_H__
