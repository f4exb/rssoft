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
 of the decoding tree

 */
 
 #ifndef _CC_TREE_GRAPHVIZ_H__
 #define _CC_TREE_GRAPHVIZ_H__
 
 #include "CC_TreeEdge.h"
 #include "CC_TreeNode.h"
 #include "CC_Encoding.h" // for print symbols
 #include <vector>
 #include <iostream>
 
 namespace ccsoft
 {
    template<typename T_IOSymbol, typename T_Register>
    class CC_TreeGraphviz
    {
    public:
        /**
         * Create a dot command sequence into an output stream
         * \param root_node Root node of the coding tree
         * \param os Output stream
         */
        static void create_dot(CC_TreeNode<T_IOSymbol, T_Register> *root_node, std::ostream& os)
        {
            std::vector<CC_TreeNode<T_IOSymbol, T_Register>*> nodes;
            std::vector<CC_TreeEdge<T_IOSymbol, T_Register>*> edges;
            
            explore_node(root_node, nodes, edges);
            print_dot(nodes, edges, os);
        }
        
    protected:
        /**
         * Explore a node in the tree to accumulate node and edges information
         * \param node Node to explore
         * \param nodes Vector of all nodes
         * \param edges Vector of all edges
         */
        static void explore_node(CC_TreeNode<T_IOSymbol, T_Register> *node,
            std::vector<CC_TreeNode<T_IOSymbol, T_Register>*>& nodes,
            std::vector<CC_TreeEdge<T_IOSymbol, T_Register>*>& edges)
        {
            nodes.push_back(node);
            const std::vector<CC_TreeEdge<T_IOSymbol, T_Register>*>& outgoing_edges = node->get_outgoing_edges();
            typename std::vector<CC_TreeEdge<T_IOSymbol, T_Register>*>::const_iterator e_it = outgoing_edges.begin();
            
            for (; e_it != outgoing_edges.end(); ++e_it)
            {
                edges.push_back(*e_it);
                explore_node((*e_it)->get_p_destination(), nodes, edges);
            }
        }
        
        /**
         * Print the dot output
         * \param nodes Vector of all nodes
         * \param edges Vector of all edges
         * \param os Output stream
         */
        static void print_dot(std::vector<CC_TreeNode<T_IOSymbol, T_Register>*>& nodes,
            std::vector<CC_TreeEdge<T_IOSymbol, T_Register>*>& edges,
            std::ostream& os)
        {
            os << "digraph G {" << std::endl;
            os << "    rankdir=LR" << std::endl << std::endl;
            
            typename std::vector<CC_TreeNode<T_IOSymbol, T_Register>*>::const_iterator n_it = nodes.begin();
            
            for (; n_it != nodes.end(); ++n_it)
            {
                os << "    n_" << (*n_it)->get_id() << " [shape=" << ((*n_it)->get_id() == 0 ? "box" : "ellipse") << ", label=\"";
                os << (*n_it)->get_id() << " " << (*n_it)->get_path_metric() << "\"";
                
                if ((*n_it)->is_on_final_path())
                {
                    os << " style=filled fillcolor=lightblue";
                }
                os  << "]" << std::endl;
            }
            
            typename std::vector<CC_TreeEdge<T_IOSymbol, T_Register>*>::const_iterator e_it = edges.begin();
            os << std::endl;
            
            for (; e_it != edges.end(); ++e_it)
            {
                os << "    n_" << (*e_it)->get_p_origin()->get_id() << " -> n_" << (*e_it)->get_p_destination()->get_id() << " [label=\"";
                print_symbol((*e_it)->get_in_symbol(), os);
                os << ":";
                print_symbol((*e_it)->get_out_symbol(),os);
                os << " " << (*e_it)->get_metric() << "\"]" << std::endl;
            }
            
            os << "}" << std::endl;
        }
    
    };
    
 } // namespace ccsoft
 
 #endif // _CC_TREE_GRAPHVIZ_H__
