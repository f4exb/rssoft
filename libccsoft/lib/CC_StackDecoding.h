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

 Convolutional soft-decision decoder based on the stack or Zigangirov-Jelinek 
 (ZJ) algorithm

 */
#ifndef __CC_STACK_DECODING_H__
#define __CC_STACK_DECODING_H__

#include "CC_Encoding.h"
#include "CCSoft_Exception.h"
#include "CC_TreeEdge.h"
#include "CC_TreeNode.h"
#include "ReliabilityMatrix.h"
#include "CC_TreeGraphviz.h"

#include <cmath>
#include <map>
#include <algorithm>


namespace ccsoft
{

float log2(float x)
{
    return log(x)/log(2.0);
}

/**
 * \brief class used for node ordering in the stack map
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


/**
 * \brief The Stack Decoding class
 * \tparam T_Register Type of the encoder internal registers
 * \tparam T_IOSymbol Type of the input and output symbols
 */
template<typename T_Register, typename T_IOSymbol>
class CC_StackDecoding
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
    CC_StackDecoding(const std::vector<unsigned int>& constraints,
            const std::vector<std::vector<T_Register> >& genpoly_representations) :
                encoding(constraints, genpoly_representations),
                root_node(0),
                use_giveup_threshold(false),
                giveup_threshold(0.0),
                use_node_limit(false),
                node_limit(0),
                codeword_score(0.0),
                cur_depth(-1),
                node_count(0),
                edge_count(0),
                tail_zeros(true)
    {}

    /**
     * Destructor. Does a final garbage collection
     */
    CC_StackDecoding()
    {
        reset();
    }
    
    /**
     * Set the give up threshold
     * \param _giveup_threshold Metric above which the process continues
     */
    void set_giveup_threshold(float _giveup_threshold)
    {
        giveup_threshold = _giveup_threshold;
        use_giveup_threshold = true;
    }
    
    /** 
     * Reset give up threshold. The process will continue until the end of the tree.
     */
    void reset_giveup_threshold()
    {
        use_giveup_threshold = false;
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
        if (root_node)
        {
            delete root_node;
            node_stack.clear();
            root_node = 0;
        }
        
        node_count = 0;
        edge_count = 0;
        codeword_score = 0.0;
        cur_depth = -1;
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
     * Get the score at the top of the stack. Valid anytime the process has started (stack not empty).
     */
    float get_stack_score() const
    {
        return node_stack.begin()->first.path_metric;
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
     * Get the stack size
     */
    unsigned int get_stack_size() const
    {
        return node_stack.size();
    }

    /**
     * Print the dot (Graphviz) file of the current decode tree to an output stream
     * \param os Output stream
     */
    void print_dot(std::ostream& os)
    {
        if (root_node)
        {
            CC_TreeGraphviz<T_IOSymbol, T_Register>::create_dot(root_node, os);
        }
    }

    /**
     * Decodes given the reliability matrix
     * \parm relmat Reference to the reliability matrix
     * \parm decoded_message Vector of symbols of retrieved message
     */
    bool decode(const ReliabilityMatrix& relmat, std::vector<T_IOSymbol>& decoded_message)
    {
        if (relmat.get_message_length() < encoding.get_m())
        {
            throw CCSoft_Exception("Reliability Matrix should have a number of columns at least equal to the code constraint");
        }

        if (relmat.get_nb_symbols_log2() != encoding.get_n())
        {
            throw CCSoft_Exception("Reliability Matrix is not compatible with code output symbol size");
        }

        reset();
        init_root(relmat); // initialize with root node
        
        // loop until we get to a terminal node or the give up condition is encountered 
        while ((node_stack.begin()->second->get_depth() < relmat.get_message_length() - 1) 
            && (!use_giveup_threshold || node_stack.begin()->second->get_path_metric() > giveup_threshold))
        {
            CC_TreeNode<T_IOSymbol, T_Register>* node = node_stack.begin()->second;
            //std::cout << std::dec << node->get_id() << ":" << node->get_depth() << ":" << node_stack.begin()->first.path_metric << std::endl;
            visit_node(node, relmat);

            if ((use_node_limit) && (node_count > node_limit))
            {
                std::cerr << "Node limit exhausted" << std::endl;
                return false;
            }
        }
        
        // Top node has the solution if we have not given up
        if (!use_giveup_threshold || node_stack.begin()->second->get_path_metric() > giveup_threshold)
        {
            //std::cout << "final: " << std::dec << node_stack.begin()->second->get_id() << ":" << node_stack.begin()->second->get_depth() << ":" << node_stack.begin()->first.path_metric << std::endl;
            back_track(node_stack.begin()->second, decoded_message, true); // back track from terminal node to retrieve decoded message
            codeword_score = node_stack.begin()->first.path_metric; // the codeword score is the path metric
            return true;
        }
        else
        {
            return false; // no solution
        }
    }

protected:
    /**
     * Initialize process at the root node
     */
    void init_root(const ReliabilityMatrix& relmat)
    {
        root_node = new CC_TreeNode<T_IOSymbol, T_Register>(node_count++, 0, 0.0, -1);
        visit_node(root_node, relmat);
    }

    /**
     * Visit a new node
     * \node Node to visit
     * \relmat Reliability matrix being used
     */
    void visit_node(CC_TreeNode<T_IOSymbol, T_Register>* node, const ReliabilityMatrix& relmat)
    {
        int forward_depth = node->get_depth() + 1;
        T_IOSymbol out_symbol;
        T_IOSymbol end_symbol;

        // return encoder to appropriate state
        if (node->get_depth() >= 0) // does not concern the root node
        {
            encoding.set_registers(node->get_registers());
        }

        if ((tail_zeros) && (forward_depth > relmat.get_message_length()-encoding.get_m()))
        {
            end_symbol = 1; // if zero tail option assume tail symbols are all zeros
        }
        else
        {
            end_symbol = (1<<encoding.get_k()); // full scan all possible input symbols
        }

        // loop through assumption for this symbol place
        for (T_IOSymbol in_symbol = 0; in_symbol < end_symbol; in_symbol++)
        {
            encoding.encode(in_symbol, out_symbol, in_symbol > 0); // step only for a new symbol place
            float edge_metric = log2(relmat(out_symbol, forward_depth));
            float forward_path_metric = edge_metric + node->get_path_metric();
            CC_TreeEdge<T_IOSymbol, T_Register> *new_edge = new CC_TreeEdge<T_IOSymbol, T_Register>(edge_count++, in_symbol, out_symbol, edge_metric, node);
            CC_TreeNode<T_IOSymbol, T_Register> *dest_node = new CC_TreeNode<T_IOSymbol, T_Register>(node_count, new_edge, forward_path_metric, forward_depth);
            dest_node->set_registers(encoding.get_registers());
            new_edge->set_p_destination(dest_node);
            node->add_outgoing_edge(new_edge); // add forward edge
            node_stack[NodeOrdering(forward_path_metric,node_count)] = dest_node;
            //std::cout << "->" << std::dec << node_count << ":" << forward_depth << " (" << (unsigned int) in_symbol << "," << (unsigned int) out_symbol << "): " << forward_path_metric << std::endl;
            node_count++;
        }
        
        cur_depth = forward_depth; // new encoder position

        if (node->get_depth() >= 0)
        {
            remove_node_from_stack(node); // remove current node from the stack unless it is the root node which is not in the stack
        }
    }
    
    /**
     * Back track from terminal node to retrieve decoded message
     */
    void back_track(CC_TreeNode<T_IOSymbol, T_Register>* node, std::vector<T_IOSymbol>& decoded_message, bool mark_nodes = false)
    {
        std::vector<T_IOSymbol> reversed_message;
        CC_TreeNode<T_IOSymbol, T_Register> *cur_node = node;
        CC_TreeEdge<T_IOSymbol, T_Register> *incoming_edge;

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
     * Removes a node from the stack map. Does a full scan but usually the nodes to be removed are on the top of the stack (i.e. beginning of the map).
     */
    void remove_node_from_stack(CC_TreeNode<T_IOSymbol, T_Register>* node)
    {
        typename std::map<NodeOrdering, CC_TreeNode<T_IOSymbol, T_Register>*, std::greater<NodeOrdering> >::iterator stack_it = node_stack.begin();
        
        for (; stack_it != node_stack.end(); ++stack_it)
        {
            if (node == stack_it->second)
            {
                node_stack.erase(stack_it);
                break;
            }
        }
    }

    CC_Encoding<T_Register, T_IOSymbol> encoding; //!< Convolutional encoding object
    CC_TreeNode<T_IOSymbol, T_Register> *root_node; //!< Root node
    std::map<NodeOrdering, CC_TreeNode<T_IOSymbol, T_Register>*, std::greater<NodeOrdering> > node_stack; //!< Ordered stack of nodes by decreasing path metric
    bool use_giveup_threshold; //!< True if a give up path metric threshold is used
    float giveup_threshold; //!< The give up path metric threshold
    bool use_node_limit; //!< Stop above number of nodes threshold
    unsigned int node_limit; //!< Number of nodes threshold
    float codeword_score; //!< Metric of the codeword found if any
    int cur_depth; //!< Current depth for the encoder
    unsigned int node_count; //!< Count of nodes in the code tree
    unsigned int edge_count; //!< Count of edges in the code tree
    bool tail_zeros; //!< True if tail of m-1 zeros in the message are assumed. This is the default option.
};

} // namespace ccsoft

#endif // __CC_STACK_DECODING_H__

