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

 Utility to output dot commands (Graphviz) for a graphical representation
 of the decoding tree. Using node+edge combination

 */
 
#ifndef _CC_TREE_GRAPHVIZ_H__
#define _CC_TREE_GRAPHVIZ_H__

#include "CC_TreeNodeEdge.h"
#include "CC_Encoding.h" // for print symbols
#include <vector>
#include <iostream>

namespace ccsoft
{
    template<typename T_IOSymbol, typename T_Register, typename T_Tag>
    class CC_TreeGraphviz
    {
    public:
        /**
         * Create a dot command sequence into an output stream
         * \param root_node Root node of the coding tree
         * \param os Output stream
         */
        static void create_dot(CC_TreeNodeEdge<T_IOSymbol, T_Register, T_Tag> *root_node, std::ostream& os)
        {
            std::vector<CC_TreeNodeEdge<T_IOSymbol, T_Register, T_Tag>*> node_edges;
            
            explore_node_edge(root_node, node_edges);
            print_dot(node_edges, os);
        }
        
    protected:
        /**
         * Explore a node+edge in the tree to accumulate node+edges information
         * \param node Node to explore
         * \param nodes Vector of all nodes
         * \param edges Vector of all edges
         */
        static void explore_node_edge(CC_TreeNodeEdge<T_IOSymbol, T_Register, T_Tag> *node_edge,
            std::vector<CC_TreeNodeEdge<T_IOSymbol, T_Register, T_Tag>*>& node_edges)
        {
            node_edges.push_back(node_edge);
            const std::vector<CC_TreeNodeEdge<T_IOSymbol, T_Register, T_Tag>*>& outgoing_node_edges = node_edge->get_outgoing_node_edges();
            typename std::vector<CC_TreeNodeEdge<T_IOSymbol, T_Register, T_Tag>*>::const_iterator ne_it = outgoing_node_edges.begin();
            
            for (; ne_it != outgoing_node_edges.end(); ++ne_it)
            {
                node_edges.push_back(*ne_it);
                explore_node_edge(*ne_it, node_edges);
            }
        }
        
        /**
         * Print the dot output
         * \param nodes Vector of all nodes
         * \param edges Vector of all edges
         * \param os Output stream
         */
        static void print_dot(std::vector<CC_TreeNodeEdge<T_IOSymbol, T_Register, T_Tag>*>& node_edges,
            std::ostream& os)
        {
            os << "digraph G {" << std::endl;
            os << "    rankdir=LR" << std::endl << std::endl;
            
            typename std::vector<CC_TreeNodeEdge<T_IOSymbol, T_Register, T_Tag>*>::const_iterator ne_it = node_edges.begin();
            
            for (; ne_it != node_edges.end(); ++ne_it)
            {
                os << "    n_" << (*ne_it)->get_id() << " [shape=" << ((*ne_it)->get_id() == 0 ? "box" : "ellipse") << ", label=\"";
                os << (*ne_it)->get_id() << " " << (*ne_it)->get_path_metric() << "\"";
                
                if ((*ne_it)->is_on_final_path())
                {
                    os << " style=filled fillcolor=lightblue";
                }
                os  << "]" << std::endl;
            }
            
            ne_it = node_edges.begin();

            for (; ne_it != node_edges.end(); ++ne_it)
            {
                os << "    n_" << (*ne_it)->get_incoming_node_edge()->get_id() << " -> n_" << (*ne_it)->get_id() << " [label=\"";
                print_symbol((*ne_it)->get_in_symbol(), os);
                os << " " << (*ne_it)->get_incoming_metric() << "\"]" << std::endl;
            }

            os << "}" << std::endl;
        }
    
    };
    
 } // namespace ccsoft
 
 #endif // _CC_TREE_GRAPHVIZ_H__
