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

#include "CC_SequentialDecoding.h"
#include "CC_Encoding.h"
#include "CCSoft_Exception.h"
#include "CC_TreeEdge.h"
#include "CC_TreeNode.h"
#include "ReliabilityMatrix.h"
#include "CC_TreeGraphviz.h"

#include <cmath>
#include <map>
#include <algorithm>
#include <iostream>


namespace ccsoft
{

/**
 * \brief The Stack Decoding class
 * \tparam T_Register Type of the encoder internal registers
 * \tparam T_IOSymbol Type of the input and output symbols
 */
template<typename T_Register, typename T_IOSymbol>
class CC_StackDecoding : public CC_SequentialDecoding<T_Register, T_IOSymbol>, public CC_SequentialDecodingInternal<T_Register, T_IOSymbol, CC_TreeEdgeTag_Empty>
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
                CC_SequentialDecoding<T_Register, T_IOSymbol>(constraints, genpoly_representations),
                CC_SequentialDecodingInternal<T_Register, T_IOSymbol, CC_TreeEdgeTag_Empty>()
    {}

    /**
     * Destructor. Does a final garbage collection
     */
    virtual ~CC_StackDecoding()
    {}

    /**
     * Reset the decoding process
     */
    void reset()
    {
        ParentInternal::reset();
        Parent::reset();
        node_stack.clear();
    }

    /**
     * Get the score at the top of the stack. Valid anytime the process has started (stack not empty).
     */
    float get_stack_score() const
    {
        return node_stack.begin()->first.path_metric;
    }

    /**
     * Get the stack size
     */
    unsigned int get_stack_size() const
    {
        return node_stack.size();
    }
    
    /**
     * Decodes given the reliability matrix
     * \param relmat Reference to the reliability matrix
     * \param decoded_message Vector of symbols of retrieved message
     */
    virtual bool decode(const ReliabilityMatrix& relmat, std::vector<T_IOSymbol>& decoded_message)
    {
        if (relmat.get_message_length() < Parent::encoding.get_m())
        {
            throw CCSoft_Exception("Reliability Matrix should have a number of columns at least equal to the code constraint");
        }

        if (relmat.get_nb_symbols_log2() != Parent::encoding.get_n())
        {
            throw CCSoft_Exception("Reliability Matrix is not compatible with code output symbol size");
        }

        reset();
        ParentInternal::init_root(relmat); // initialize the root node
        Parent::node_count++;
        visit_node_forward(ParentInternal::root_node, relmat); // visit the root node

        // loop until we get to a terminal node or the metric limit is encountered hence the stack is empty
        while ((node_stack.size() > 0)
            && (node_stack.begin()->second->get_depth() < relmat.get_message_length() - 1))
        {
            StackNode* node = node_stack.begin()->second;
            //std::cout << std::dec << node->get_id() << ":" << node->get_depth() << ":" << node_stack.begin()->first.path_metric << std::endl;
            visit_node_forward(node, relmat);

            if ((Parent::use_node_limit) && (Parent::node_count > Parent::node_limit))
            {
                std::cerr << "Node limit exhausted" << std::endl;
                return false;
            }
        }

        // Top node has the solution if we have not given up
        if (!Parent::use_metric_limit || node_stack.size() != 0)
        {
            //std::cout << "final: " << std::dec << node_stack.begin()->second->get_id() << ":" << node_stack.begin()->second->get_depth() << ":" << node_stack.begin()->first.path_metric << std::endl;
            ParentInternal::back_track(node_stack.begin()->second, decoded_message, true); // back track from terminal node to retrieve decoded message
            Parent::codeword_score = node_stack.begin()->first.path_metric; // the codeword score is the path metric
            return true;
        }
        else
        {
            std::cerr << "Metric limit encountered" << std::endl;
            return false; // no solution
        }
    }
    
    /**
     * Print stats to an output stream
     * \param os Output stream
     * \param success True if decoding was successful
     */
    virtual void print_stats(std::ostream& os, bool success)
    {
        std::cout << "score = " << Parent::get_score()
                << " stack_score = " << get_stack_score()
                << " #nodes = " << Parent::get_nb_nodes()
                << " stack_size = " << get_stack_size()
                << " max depth = " << Parent::get_max_depth() << std::endl;
        std::cout << "_RES " << (success ? 1 : 0) << ","
                << Parent::get_score() << ","
                << get_stack_score() << ","
                << Parent::get_nb_nodes() << ","
                << get_stack_size() << ","
                << Parent::get_max_depth() << std::endl;
    }
    
    /**
     * Print the dot (Graphviz) file of the current decode tree to an output stream
     * \param os Output stream
     */
    virtual void print_dot(std::ostream& os)
    {
        ParentInternal::print_dot_internal(os);
    }

protected:
    typedef CC_SequentialDecoding<T_Register, T_IOSymbol> Parent;                                       //!< Parent class this class inherits from
    typedef CC_SequentialDecodingInternal<T_Register, T_IOSymbol, CC_TreeEdgeTag_Empty> ParentInternal; //!< Parent class this class inherits from
    typedef CC_TreeNode<T_IOSymbol, T_Register, CC_TreeEdgeTag_Empty> StackNode; //!< Class of code tree nodes in the stack algorithm
    typedef CC_TreeEdge<T_IOSymbol, T_Register, CC_TreeEdgeTag_Empty> StackEdge; //!< Class of code tree edges in the stack algorithm

    /**
     * Visit a new node
     * \node Node to visit
     * \relmat Reliability matrix being used
     */
    virtual void visit_node_forward(CC_TreeNode<T_IOSymbol, T_Register, CC_TreeEdgeTag_Empty>* node, const ReliabilityMatrix& relmat)
    {
        int forward_depth = node->get_depth() + 1;
        T_IOSymbol out_symbol;
        T_IOSymbol end_symbol;

        // return encoder to appropriate state
        if (node->get_depth() >= 0) // does not concern the root node
        {
            Parent::encoding.set_registers(node->get_registers());
        }

        if ((Parent::tail_zeros) && (forward_depth > relmat.get_message_length()-Parent::encoding.get_m()))
        {
            end_symbol = 1; // if zero tail option assume tail symbols are all zeros
        }
        else
        {
            end_symbol = (1<<Parent::encoding.get_k()); // full scan all possible input symbols
        }

        // loop through assumption for this symbol place
        for (T_IOSymbol in_symbol = 0; in_symbol < end_symbol; in_symbol++)
        {
            Parent::encoding.encode(in_symbol, out_symbol, in_symbol > 0); // step only for a new symbol place
            float edge_metric = log2(relmat(out_symbol, forward_depth)) - Parent::edge_bias;
            
            float forward_path_metric = edge_metric + node->get_path_metric();
            if ((!Parent::use_metric_limit) || (forward_path_metric > Parent::metric_limit))
            {
                StackEdge *new_edge = new StackEdge(Parent::edge_count++, in_symbol, out_symbol, edge_metric, node);
                StackNode *dest_node = new StackNode(Parent::node_count, new_edge, forward_path_metric, forward_depth);
                dest_node->set_registers(Parent::encoding.get_registers());
                new_edge->set_p_destination(dest_node);
                node->add_outgoing_edge(new_edge); // add forward edge
                node_stack[NodeOrdering(forward_path_metric, Parent::node_count)] = dest_node;
                //std::cout << "->" << std::dec << node_count << ":" << forward_depth << " (" << (unsigned int) in_symbol << "," << (unsigned int) out_symbol << "): " << forward_path_metric << std::endl;
                Parent::node_count++;
            }
        }

        Parent::cur_depth = forward_depth; // new encoder position

        if (Parent::cur_depth > Parent::max_depth)
        {
        	Parent::max_depth = Parent::cur_depth;
        }

        if (node->get_depth() >= 0)
        {
            remove_node_from_stack(node); // remove current node from the stack unless it is the root node which is not in the stack
        }
    }

    /**
     * Removes a node from the stack map. Does a full scan but usually the nodes to be removed are on the top of the stack (i.e. beginning of the map).
     */
    void remove_node_from_stack(StackNode* node)
    {
        typename std::map<NodeOrdering, StackNode*, std::greater<NodeOrdering> >::iterator stack_it = node_stack.begin();

        for (; stack_it != node_stack.end(); ++stack_it)
        {
            if (node == stack_it->second)
            {
                node_stack.erase(stack_it);
                break;
            }
        }
    }
    
    std::map<NodeOrdering, StackNode*, std::greater<NodeOrdering> > node_stack; //!< Ordered stack of nodes by decreasing path metric
};

} // namespace ccsoft

#endif // __CC_STACK_DECODING_H__

