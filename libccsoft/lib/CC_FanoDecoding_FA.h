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
 web.ntpu.edu.tw/~yshan/book_chapter.pdf

 Uses a node+edge combo code tree representation‎

 */
#ifndef __CC_FANO_DECODING_FA_H__
#define __CC_FANO_DECODING_FA_H__

#include "CC_SequentialDecoding_FA.h"
#include "CC_SequentialDecodingInternal_FA.h"
#include "CC_Encoding_FA.h"
#include "CCSoft_Exception.h"
#include "CC_TreeNodeEdge_FA.h"
#include "CC_ReliabilityMatrix.h"
#include "Debug.h"

#include <cmath>
#include <time.h>
#include <algorithm>
#include <iostream>
#include <iomanip>

namespace ccsoft
{

/**
 * \brief The Fano like Decoding class. Tag is a boolean used as the traversed back indicator.
 * This version uses fixed arrays to store registers and forward node+edges pointers.
 * N_k template parameter gives the size of the input symbol (k parameter) and therefore the number of registers.
 * There are (1<<N_k) forward node+edges.
 * \tparam T_Register Type of the encoder internal registers
 * \tparam T_IOSymbol Type of the input and output symbols
 * \tparam N_k Size of an input symbol in bits (k parameter)
 */
template<typename T_Register, typename T_IOSymbol, unsigned int N_k>
class CC_FanoDecoding_FA : public CC_SequentialDecoding_FA<T_Register, T_IOSymbol, N_k> ,public CC_SequentialDecodingInternal_FA<T_Register, T_IOSymbol, bool, N_k>
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
	CC_FanoDecoding_FA(const std::vector<unsigned int>& constraints,
            const std::vector<std::vector<T_Register> >& genpoly_representations,
            float _init_threshold,
            float _delta_threshold,
            unsigned int _tree_cache_size = 0,
            float _delta_init_threshold = 0.0) :
                CC_SequentialDecoding_FA<T_Register, T_IOSymbol, N_k>(constraints, genpoly_representations),
                CC_SequentialDecodingInternal_FA<T_Register, T_IOSymbol, bool, N_k>(),
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
    virtual ~CC_FanoDecoding_FA()
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
        FanoNodeEdge *node_edge_current, *node_edge_successor;

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
        node_edge_current = ParentInternal::root_node;
        nb_moves = 0;

#ifdef _DEBUG
        timespec time1, time2;
        int time_option = CLOCK_REALTIME;
        clock_gettime(time_option, &time1);
#endif

        visit_node_forward(node_edge_current, relmat);

        while (continue_process(node_edge_current, relmat))
        {
            //std::cout << "T=" << cur_threshold << " depth=" << node_current->get_depth() << " node #" << node_current->get_id() << " Mc=" << node_current->get_path_metric() << std::endl;
            DEBUG_OUT(Parent::verbosity > 1, "T=" << cur_threshold << " depth=" << node_edge_current->get_depth() << " node #" << node_edge_current->get_id() << " Mc=" << node_edge_current->get_path_metric() << std::endl);

            if (node_edge_current->get_depth() > Parent::max_depth)
            {
            	Parent::max_depth = node_edge_current->get_depth();
            }

            if (node_edge_current == ParentInternal::root_node)
            {
                root_threshold = cur_threshold;
            }

            nb_moves++;
            const std::array<FanoNodeEdge*, (1<<N_k)>& outgoing_node_edges = node_edge_current->get_outgoing_node_edges();
            typename std::array<FanoNodeEdge*, (1<<N_k)>::const_iterator ne_it = outgoing_node_edges.begin();
            std::vector<FanoNodeEdge*> child_node_edges;

            for (; ne_it != outgoing_node_edges.end(); ++ne_it)
            {
                if ((*ne_it) && !((*ne_it)->get_tag())) // not traversed back
                {
                    child_node_edges.push_back(*ne_it);
                }
            }

            if (child_node_edges.size() == 0) // exhausted forward paths
            {
                DEBUG_OUT(Parent::verbosity > 2, "exhaustion of forward paths at node #" << node_edge_current->get_id() << std::endl);
                node_edge_current = move_back_from_node_or_loosen_threshold(node_edge_current);
                continue;
            }

            std::sort(child_node_edges.begin(), child_node_edges.end(), node_edge_pointer_ordering<FanoNodeEdge>);
            node_edge_successor = *child_node_edges.begin(); // best successor
            DEBUG_OUT(Parent::verbosity > 2, "best successor node #" << node_edge_successor->get_id() << " Ms=" << node_edge_successor->get_path_metric() << std::endl);

            if (node_edge_successor->get_path_metric() >= cur_threshold) // Ms >= T
            {
                // move forward:
                DEBUG_OUT(Parent::verbosity > 2, "forward" << std::endl);
                FanoNodeEdge *node_predecessor = node_edge_current;
                node_edge_current = node_edge_successor;

                // termination with solution
                if (node_edge_current->get_depth() == relmat.get_message_length() - 1)
                {
                	Parent::codeword_score = node_edge_current->get_path_metric();
                    ParentInternal::back_track(node_edge_current, decoded_message, true); // back track from terminal node to retrieve decoded message
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
                    int nb_delta = int((node_edge_current->get_path_metric() - init_threshold) / delta_threshold);

                    if (nb_delta < 0)
                    {
                    	cur_threshold = ((nb_delta - 1) * delta_threshold) + init_threshold;
                    }
                    else
                    {
                    	cur_threshold = (nb_delta * delta_threshold) + init_threshold;
                    }

                    DEBUG_OUT(Parent::verbosity > 2, "tightening " << node_edge_current->get_path_metric() << " -> " << cur_threshold << std::endl);
                }

                // create children nodes from the new current node
                visit_node_forward(node_edge_current, relmat);
            }
            else
            {
                node_edge_current = move_back_from_node_or_loosen_threshold(node_edge_current);
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
    typedef CC_SequentialDecoding_FA<T_Register, T_IOSymbol, N_k> Parent; //!< Parent class this class inherits from
    typedef CC_SequentialDecodingInternal_FA<T_Register, T_IOSymbol, bool, N_k> ParentInternal; //!< Parent class this class inherits from
    typedef CC_TreeNodeEdge_FA<T_IOSymbol, T_Register, bool, N_k> FanoNodeEdge;   //!< Class of code tree nodes in the Fano algorithm

    /**
     * Visit a new node
     * \parm node Node to visit
     * \parm relmat Reliability matrix being used
     */
    virtual void visit_node_forward(FanoNodeEdge* node_edge, const CC_ReliabilityMatrix& relmat)
    {
        unsigned int n = Parent::encoding.get_n();
        int forward_depth = node_edge->get_depth() + 1;
        T_IOSymbol out_symbol;
        T_IOSymbol end_symbol;

        // return encoder to appropriate state
        if (node_edge->get_depth() >= 0) // does not concern the root node_edge
        {
            Parent::encoding.set_registers(node_edge->get_registers());
        }

        if ((Parent::tail_zeros) && (forward_depth > relmat.get_message_length()-Parent::encoding.get_m()))
        {
            end_symbol = 1; // if zero tail option assume tail symbols are all zeros
        }
        else
        {
            end_symbol = (1<<Parent::encoding.get_k()); // full scan all possible input symbols
        }

        if (!node_edge->valid_outgoing_node_edges(end_symbol)) // edges are not cached
        {
            if ((tree_cache_size > 0) && (effective_node_count >= tree_cache_size)) // if tree cache is used and cache limit reached
            {
                purge_tree_cache(node_edge); // purge before allocating new nodes
            }

            // loop through assumption for this symbol place and create child nodes
            for (T_IOSymbol in_symbol = 0; in_symbol < end_symbol; in_symbol++)
            {
                Parent::encoding.encode(in_symbol, out_symbol, in_symbol > 0); // step only for a new symbol place
                float edge_metric = ParentInternal::log2(relmat(out_symbol, forward_depth)) - Parent::edge_bias;
                float forward_path_metric = edge_metric + node_edge->get_path_metric();
                FanoNodeEdge *next_node_edge = new FanoNodeEdge(Parent::node_count++, node_edge, in_symbol, edge_metric, forward_path_metric, forward_depth);
                next_node_edge->get_tag() = false; // Init traversed back indicator
                next_node_edge->set_registers(Parent::encoding.get_registers());
                node_edge->set_outgoing_node_edge(next_node_edge, in_symbol); // add forward edge
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
    FanoNodeEdge *move_back_from_node_or_loosen_threshold(FanoNodeEdge *node_edge_current)
    {
        if (node_edge_current == ParentInternal::root_node) // at root node there are no other options than loosening threshold
        {
            cur_threshold -= delta_threshold;
            DEBUG_OUT(Parent::verbosity > 2, "loosening " << node_edge_current->get_path_metric() << " -> " << cur_threshold << std::endl);
        }
        else
        {
            FanoNodeEdge *node_edge_predecessor = node_edge_current->get_incoming_node_edge();

            if (node_edge_predecessor->get_path_metric() >= cur_threshold) // move backward
            {
                DEBUG_OUT(Parent::verbosity > 2, std::cout << "backward" << std::endl);

                if (tree_cache_size == 0) // tree cache is not used
                {
                    // delete all successor edges and nodes
                    std::array<FanoNodeEdge*, (1<<N_k)>& outgoing_node_edges = node_edge_current->get_outgoing_node_edges();
                    typename std::array<FanoNodeEdge*, (1<<N_k)>::iterator ne_it = outgoing_node_edges.begin();

                    for (;ne_it != outgoing_node_edges.end(); ++ne_it)
                    {
                        delete *ne_it;
                    }

                    effective_node_count -= outgoing_node_edges.size();
                    outgoing_node_edges.fill(0);
                }

                // mark incoming edge as traversed back
                if (node_edge_predecessor != ParentInternal::root_node)
                {
                    node_edge_current->get_tag() = true;
                }

                // move back: change node address to previous node address
                node_edge_current = node_edge_predecessor;
            }
            else // loosen threshold
            {
                cur_threshold -= delta_threshold;
                DEBUG_OUT(Parent::verbosity > 2, "loosening " << node_edge_current->get_path_metric() << " -> " << cur_threshold << std::endl);
            }
        }

        return node_edge_current;
    }

    /**
     * Check if process can continue
     */
    bool continue_process(FanoNodeEdge *node_edge_current, const CC_ReliabilityMatrix& relmat)
    {
        if ((node_edge_current == ParentInternal::root_node) && (nb_moves > 0) && (cur_threshold == root_threshold))
        {
            const std::array<FanoNodeEdge*, (1<<N_k)>& outgoing_node_edges = node_edge_current->get_outgoing_node_edges();
            typename std::array<FanoNodeEdge*, (1<<N_k)>::const_iterator ne_it = outgoing_node_edges.begin();
            bool children_open = true;

            for (; ne_it != outgoing_node_edges.end(); ++ne_it)
            {
                if ((*ne_it)->get_tag()) // traversed back
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
                    ParentInternal::root_node->delete_outgoing_node_edges(); // effectively resets the root node without destroying it
                    Parent::node_count = 1;
                    effective_node_count = 1;
                    nb_moves = 0;
                    visit_node_forward(node_edge_current, relmat); // visit root node again
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
     * Purge tree cache from a node. Keep only the current path to the node and path nodes immediate successors
     * \param node The node to purge from
     */
    void purge_tree_cache(FanoNodeEdge *node_edge)
    {
        bool node_terminal = true;
        unsigned int remaining_nodes = 0;

        while (node_edge != ParentInternal::root_node)
        {
            FanoNodeEdge *node_edge_predecessor = node_edge->get_incoming_node_edge();
            std::array<FanoNodeEdge*, (1<<N_k)>& outgoing_node_edges = node_edge_predecessor->get_outgoing_node_edges();
            typename std::array<FanoNodeEdge*, (1<<N_k)>::iterator ne_it = outgoing_node_edges.begin();

            for (;ne_it != outgoing_node_edges.end(); ++ne_it)
            {
                FanoNodeEdge *node_edge_sibling = (*ne_it);

                if (node_terminal || (node_edge_sibling != node_edge))
                {
                    (*ne_it)->delete_outgoing_node_edges();
                }

                remaining_nodes++;
            }

            node_edge = node_edge_predecessor;
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

#endif // __CC_FANO_DECODING_FA_H__
