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

 Convolutional soft-decision decoder based on the Fano sequential algorithm as described in
 Sequential Decoding of Convolutional Codes by Yunghsiang S. Han and Po-Ning Chen (algorithm p.26)
 web.ntpu.edu.tw/~yshan/book_chapter.pdf‎

 */
#ifndef __CC_FANO_DECODING_H__
#define __CC_FANO_DECODING_H__

#include "CC_SequentialDecoding.h"
#include "CC_Encoding.h"
#include "CCSoft_Exception.h"
#include "CC_TreeEdge.h"
#include "CC_TreeNode.h"
#include "CC_ReliabilityMatrix.h"
#include "CC_TreeGraphviz.h"
#include "Debug.h"

#include <cmath>
#include <time.h>
#include <algorithm>
#include <iostream>
#include <iomanip>

namespace ccsoft
{

/**
 * \brief The Fano like Decoding class. Edge tag is a boolean used as the traversed back indicator.
 * \tparam T_Register Type of the encoder internal registers
 * \tparam T_IOSymbol Type of the input and output symbols
 */
template<typename T_Register, typename T_IOSymbol>
class CC_FanoDecoding : public CC_SequentialDecoding<T_Register, T_IOSymbol> ,public CC_SequentialDecodingInternal<T_Register, T_IOSymbol, bool>
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
     * \param _init_threshold Initial path metric threshold
     * \param _delta_threshold Delta of path metric that is applied when lowering threshold
     * \param _tree_cache_size Tree cache maximum size in number of nodes (0 if not used)
     * \param _delta_init_threshold: Delta of path metric that is applied when restarting with a lower initial threshold (0 if not used)
     */
	CC_FanoDecoding(const std::vector<unsigned int>& constraints,
            const std::vector<std::vector<T_Register> >& genpoly_representations,
            float _init_threshold,
            float _delta_threshold,
            unsigned int _tree_cache_size = 0,
            float _delta_init_threshold = 0.0) :
                CC_SequentialDecoding<T_Register, T_IOSymbol>(constraints, genpoly_representations),
                CC_SequentialDecodingInternal<T_Register, T_IOSymbol, bool>(),
                init_threshold(_init_threshold),
                cur_threshold(_init_threshold),
                root_threshold(_init_threshold),
                delta_threshold(_delta_threshold),
                solution_found(false),
                effective_node_count(0),
                nb_moves(0),
                tree_cache_size(_tree_cache_size),
                unloop(_delta_init_threshold < 0.0),
                delta_init_threshold(_delta_init_threshold)
    {}

    /**
     * Destructor. Does a final garbage collection
     */
    virtual ~CC_FanoDecoding()
    {}

    /**
     * Set the tree cache size
     * \param _tree_cache_size Maximum number of nodes to be cached
     */
    void set_tree_cache_size(unsigned int _tree_cache_size)
    {
        tree_cache_size = _tree_cache_size;
    }
    
    /**
     * Reset the decoding process
     */
    void reset()
    {
        ParentInternal::reset();
        Parent::reset();
        cur_threshold = init_threshold;
        solution_found = false;
        effective_node_count = 0;
    }

    /**
     * Decodes given the reliability matrix. Algorithm reproduced from Sequential Decoding of Convolutional Codes
     * by Yunghsiang S. Han and Po-Ning Chen p.26
     * \parm relmat Reference to the reliability matrix
     * \parm decoded_message Vector of symbols of retrieved message
     */
    virtual bool decode(const CC_ReliabilityMatrix& relmat, std::vector<T_IOSymbol>& decoded_message)
    {
        FanoNode *node_current, *node_successor;

        if (relmat.get_message_length() < Parent::encoding.get_m())
        {
            throw CCSoft_Exception("Reliability Matrix should have a number of columns at least equal to the code constraint");
        }

        if (relmat.get_nb_symbols_log2() != Parent::encoding.get_n())
        {
            throw CCSoft_Exception("Reliability Matrix is not compatible with code output symbol size");
        }

        reset();
        ParentInternal::init_root(); // initialize root node
        Parent::node_count++;
        effective_node_count++;
        node_current = ParentInternal::root_node;
        nb_moves = 0;

#ifdef _DEBUG
        timespec time1, time2;
        int time_option = CLOCK_REALTIME;
        clock_gettime(time_option, &time1);
#endif

        visit_node_forward(node_current, relmat);

        while (continue_process(node_current, relmat))
        {
            //std::cout << "T=" << cur_threshold << " depth=" << node_current->get_depth() << " node #" << node_current->get_id() << " Mc=" << node_current->get_path_metric() << std::endl;
            DEBUG_OUT(Parent::verbosity > 1, "T=" << cur_threshold << " depth=" << node_current->get_depth() << " node #" << node_current->get_id() << " Mc=" << node_current->get_path_metric() << std::endl);

            if (node_current->get_depth() > Parent::max_depth)
            {
            	Parent::max_depth = node_current->get_depth();
            }

            if (node_current == ParentInternal::root_node)
            {
                root_threshold = cur_threshold;
            }

            nb_moves++;
            const std::vector<FanoEdge*>& outgoing_edges = node_current->get_outgoing_edges();
            typename std::vector<FanoEdge*>::const_iterator e_it = outgoing_edges.begin();
            std::vector<FanoNode*> child_nodes;

            for (; e_it != outgoing_edges.end(); ++e_it)
            {
                if (!((*e_it)->get_edge_tag())) // not traversed back
                {
                    child_nodes.push_back((*e_it)->get_p_destination());
                }
            }

            if (child_nodes.size() == 0) // exhausted forward paths
            {
                DEBUG_OUT(Parent::verbosity > 2, "exhaustion of forward paths at node #" << node_current->get_id() << std::endl);
                node_current = move_back_from_node_or_loosen_threshold(node_current);
                continue;
            }

            std::sort(child_nodes.begin(), child_nodes.end(), node_pointer_ordering<FanoNode>);
            node_successor = *child_nodes.begin(); // best successor
            DEBUG_OUT(Parent::verbosity > 2, "best successor node #" << node_successor->get_id() << " Ms=" << node_successor->get_path_metric() << std::endl);

            if (node_successor->get_path_metric() >= cur_threshold) // Ms >= T
            {
                // move forward:
                DEBUG_OUT(Parent::verbosity > 2, "forward" << std::endl);
                FanoNode *node_predecessor = node_current;
                node_current = node_successor;

                // termination with solution
                if (node_current->get_depth() == relmat.get_message_length() - 1)
                {
                	Parent::codeword_score = node_current->get_path_metric();
                    ParentInternal::back_track(node_current, decoded_message, true); // back track from terminal node to retrieve decoded message
                    solution_found = true;
                    Parent::max_depth++;
#ifdef _DEBUG
                    clock_gettime(time_option, &time2);
                    DEBUG_OUT(Parent::verbosity > 0, std::cout << "Decoding time: " << std::setw(12) << std::setprecision(9) << debug_get_time_difference(time2,time1) << " s" << std::endl);
#endif
                    return true;
                }

                // threshold tightening for the new current node
                if (node_predecessor->get_path_metric() < cur_threshold + delta_threshold)
                {
                    int nb_delta = int((node_current->get_path_metric() - init_threshold) / delta_threshold);

                    if (nb_delta < 0)
                    {
                    	cur_threshold = ((nb_delta - 1) * delta_threshold) + init_threshold;
                    }
                    else
                    {
                    	cur_threshold = (nb_delta * delta_threshold) + init_threshold;
                    }

                    DEBUG_OUT(Parent::verbosity > 2, "tightening " << node_current->get_path_metric() << " -> " << cur_threshold << std::endl);
                }

                // create children nodes from the new current node
                visit_node_forward(node_current, relmat);
            }
            else
            {
                node_current = move_back_from_node_or_loosen_threshold(node_current);
            }
        }

        return false;
    }

    /**
     * Print stats to an output stream
     * \param os Output stream
     * \param success True if decoding was successful
     */
    virtual void print_stats(std::ostream& os, bool success)
    {
        std::cout << "score = " << Parent::get_score()
                << " cur.threshold = " << cur_threshold
                << " nodes = " << Parent::get_nb_nodes()
                << " eff.nodes = " << effective_node_count
                << " moves = " << nb_moves
                << " max depth = " << Parent::get_max_depth() << std::endl;
        std::cout << "_RES " << (success ? 1 : 0) << ","
                << Parent::get_score() << ","
                << cur_threshold << ","
                << Parent::get_nb_nodes() << ","
                << effective_node_count << ","
                << nb_moves << ","
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
    typedef CC_SequentialDecoding<T_Register, T_IOSymbol> Parent; //!< Parent class this class inherits from
    typedef CC_SequentialDecodingInternal<T_Register, T_IOSymbol, bool> ParentInternal; //!< Parent class this class inherits from
    typedef CC_TreeNode<T_IOSymbol, T_Register, bool> FanoNode;   //!< Class of code tree nodes in the Fano algorithm
    typedef CC_TreeEdge<T_IOSymbol, T_Register, bool> FanoEdge;   //!< Class of code tree edges in the Fano algorithm

    /**
     * Visit a new node
     * \parm node Node to visit
     * \parm relmat Reliability matrix being used
     */
    virtual void visit_node_forward(CC_TreeNode<T_IOSymbol, T_Register, bool>* node, const CC_ReliabilityMatrix& relmat)
    {
        unsigned int n = Parent::encoding.get_n();
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

        if (node->get_outgoing_edges().size() == 0) // edges are not cached
        {
            if ((tree_cache_size > 0) && (effective_node_count >= tree_cache_size)) // if tree cache is used and cache limit reached
            {
                purge_tree_cache(node); // purge before allocating new nodes
            }

            // loop through assumption for this symbol place and create child nodes
            for (T_IOSymbol in_symbol = 0; in_symbol < end_symbol; in_symbol++)
            {
                Parent::encoding.encode(in_symbol, out_symbol, in_symbol > 0); // step only for a new symbol place
                float edge_metric = log2(relmat(out_symbol, forward_depth)) - Parent::edge_bias;
                float forward_path_metric = edge_metric + node->get_path_metric();
                FanoEdge *new_edge = new FanoEdge(Parent::edge_count++, in_symbol, out_symbol, edge_metric, node);
                new_edge->get_edge_tag() = false; // Init traversed back indicator
                FanoNode *dest_node = new FanoNode(Parent::node_count++, new_edge, forward_path_metric, forward_depth);
                dest_node->set_registers(Parent::encoding.get_registers());
                new_edge->set_p_destination(dest_node);
                node->add_outgoing_edge(new_edge); // add forward edge
                effective_node_count++;
            }
        }
    }

    /**
     * Chooses between moving back from the node or loosen threshold
     * Before moving back it deletes all successors of the node (edges and nodes) and
     * it marks the incoming edge as traversed back
     * \param node_current Node to move backfrom or where to loosen threshold
     */
    FanoNode *move_back_from_node_or_loosen_threshold(FanoNode *node_current)
    {
        if (node_current == ParentInternal::root_node) // at root node there are no other options than loosening threshold
        {
            cur_threshold -= delta_threshold;
            DEBUG_OUT(Parent::verbosity > 2, "loosening " << node_current->get_path_metric() << " -> " << cur_threshold << std::endl);
        }
        else
        {
            FanoNode *node_predecessor = node_current->get_incoming_edge()->get_p_origin();

            if (node_predecessor->get_path_metric() >= cur_threshold) // move backward
            {
                DEBUG_OUT(Parent::verbosity > 2, std::cout << "backward" << std::endl);

                if (tree_cache_size == 0) // tree cache is not used
                {
                    // delete all successor edges and nodes
                    std::vector<FanoEdge*>& outgoing_edges = node_current->get_outgoing_edges();
                    typename std::vector<FanoEdge*>::iterator e_it = outgoing_edges.begin();

                    for (;e_it != outgoing_edges.end(); ++e_it)
                    {
                        delete *e_it;
                    }

                    effective_node_count -= outgoing_edges.size();
                    outgoing_edges.clear();
                }

                // mark incoming edge as traversed back
                if (node_predecessor != ParentInternal::root_node)
                {
                    node_current->get_incoming_edge()->get_edge_tag() = true;
                }

                // move back: change node address to previous node address
                node_current = node_predecessor;
            }
            else // loosen threshold
            {
                cur_threshold -= delta_threshold;
                DEBUG_OUT(Parent::verbosity > 2, "loosening " << node_current->get_path_metric() << " -> " << cur_threshold << std::endl);
            }
        }

        return node_current;
    }

    /**
     * Check if process can continue
     */
    bool continue_process(FanoNode *node_current, const CC_ReliabilityMatrix& relmat) 
    {
        if ((node_current == ParentInternal::root_node) && (nb_moves > 0) && (cur_threshold == root_threshold))
        {
            const std::vector<FanoEdge*>& outgoing_edges = node_current->get_outgoing_edges();
            typename std::vector<FanoEdge*>::const_iterator e_it = outgoing_edges.begin();
            bool children_open = true;

            for (; e_it != outgoing_edges.end(); ++e_it)
            {
                if ((*e_it)->get_edge_tag()) // traversed back
                {
                    children_open = false;
                    break;
                }
            }

            if (children_open)
            {
                if (unloop && ((Parent::use_metric_limit) && (init_threshold > Parent::metric_limit)))
                {
                    init_threshold += delta_init_threshold; // lower initial threshold and start all over again (delta if used is negative)
                    Parent::reset();                        // reset but do not delete root node
                    cur_threshold = init_threshold;
                    solution_found = false;
                    ParentInternal::root_node->delete_outgoing_edges(); // effectively resets the root node without destroying it
                    Parent::node_count = 1;
                    effective_node_count = 1;
                    nb_moves = 0;
                    visit_node_forward(node_current, relmat); // visit root node again
                    std::cerr << "Loop condition detected, restart with init threshold = " << init_threshold << std::endl;
                    return true;
                }
                else
                {
                    std::cerr << "Loop condition detected, aborting" << std::endl;
                    return false;
                }
            }
        }

        if ((Parent::use_metric_limit) && (cur_threshold < Parent::metric_limit))
        {
            std::cerr << "Metric limit encountered" << std::endl;
            return false;
        }

        if ((Parent::use_node_limit) && (Parent::node_count > Parent::node_limit))
        {
            std::cerr << "Node limit exhausted" << std::endl;
            return false;
        }

        return true;
    }

    /**
     * Purge tree cache from a node. Keep only the current path and path nodes immediate successors
     * \param node The node to purge from
     */
    void purge_tree_cache(FanoNode *node)
    {
        bool node_terminal = true;
        unsigned int remaining_nodes = 0;

        while (node != ParentInternal::root_node)
        {
            FanoNode *node_predecessor = node->get_incoming_edge()->get_p_origin();
            std::vector<FanoEdge*>& outgoing_edges = node_predecessor->get_outgoing_edges();
            typename std::vector<FanoEdge*>::iterator e_it = outgoing_edges.begin();

            for (;e_it != outgoing_edges.end(); ++e_it)
            {
                FanoNode *node_sibling = (*e_it)->get_p_destination();

                if (node_terminal || (node_sibling != node))
                {
                    (*e_it)->get_p_destination()->delete_outgoing_edges();
                }

                remaining_nodes++;
            }

            node = node_predecessor;
            node_terminal = false;
        }

        remaining_nodes++; // +1 for root node
        effective_node_count = remaining_nodes;
        DEBUG_OUT(Parent::verbosity > 1, "purged tree cache, nb of remaining nodes = " << remaining_nodes << std::endl);
    }

    float init_threshold;              //!< Initial path metric threshold
    float cur_threshold;               //!< Current path metric threshold
    float delta_threshold;             //!< Delta of path metric that is applied when lowering threshold
    bool solution_found;               //!< Set to true when eligible terminal node is found
    unsigned int effective_node_count; //!< Count of nodes effectively present in the system
    unsigned int nb_moves;             //!< Number of moves i.e. number of iterations in the main loop
    float root_threshold;              //!< Latest threshold at root node
    unsigned int tree_cache_size;      //!< Tree cache size in maximum number of nodes in cache (0 = tree is not cached)
    bool unloop;                       //!< If true when a loop condition is detected attempt to restart with a lower threshold
    float delta_init_threshold;        //!< Delta of path metric that is applied when restarting with a lower initial threshold 
};



} // namespace ccsoft

#endif // __CC_FANO_DECODING_H__
