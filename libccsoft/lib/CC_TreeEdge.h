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

 Edge in the convolutional code tree

 */
#ifndef __CC_TREE_EDGE_H__
#define __CC_TREE_EDGE_H__

namespace ccsoft
{

class CC_TreeEdgeTag_Empty
{
};

template<typename T_IOSymbol, typename T_Register, typename T_EdgeTag>
class CC_TreeNode;

/**
 * \brief An edge of the code tree
 * \tparam T_IOSymbol Type of the input and output symbols
 * \tparam T_Register Type of the encoder internal registers
 * \tparam T_Register Type of the edge tag
 */
template<typename T_IOSymbol, typename T_Register, typename T_EdgeTag>
class CC_TreeEdge
{
public:
    /**
     * Constructor
     * \param _id Unique ID of the edge
     * \param _in_symbol Input symbol corresponding to the edge
     * \param _out_symbol Output symbol corresponding to the edge
     * \param _metric Metric of the edge
     * \param _p_origin Pointer to the edge origin node
     */
    CC_TreeEdge(unsigned int _id,
            const T_IOSymbol& _in_symbol,
            const T_IOSymbol& _out_symbol,
            float _metric,
            CC_TreeNode<T_IOSymbol, T_Register, T_EdgeTag> *_p_origin) :
                id(_id),
                in_symbol(_in_symbol),
                out_symbol(_out_symbol),
                metric(_metric),
                p_origin(_p_origin),
                p_destination(0)
    {}

    /**
     * Destructor. Destroys the destination node
     */
    ~CC_TreeEdge()
    {
        if (p_destination)
        {
            delete p_destination;
            p_destination = 0;
        }
    }
    
    /** 
     * Edge node destination setter
     */
    void set_p_destination(CC_TreeNode<T_IOSymbol, T_Register, T_EdgeTag> *_p_destination)
    {
        p_destination = _p_destination;
    }

    /**
     * Input symbol getter
     */
    const T_IOSymbol& get_in_symbol() const
    {
        return in_symbol;
    }

    /**
     * Output symbol getter
     */
    const T_IOSymbol& get_out_symbol() const
    {
        return out_symbol;
    }

    /**
     * Metric getter
     */
    float get_metric() const
    {
        return metric;
    }

    /**
     * Origin pointer getter
     */
    CC_TreeNode<T_IOSymbol, T_Register, T_EdgeTag> *get_p_origin()
    {
        return p_origin;
    }

    /**
     * Destination pointer getter
     */
    CC_TreeNode<T_IOSymbol, T_Register, T_EdgeTag> *get_p_destination()
    {
        return p_destination;
    }
    
    /**
     * R/O reference to edge tag
     */
    const T_EdgeTag& get_edge_tag() const
    {
        return edge_tag;
    }

    /**
     * R/W reference to edge tag
     */
    T_EdgeTag& get_edge_tag()
    {
        return edge_tag;
    }

protected:
    unsigned int id; //!< Edge's unique ID
    T_IOSymbol in_symbol; //!< Input symbol corresponding to the edge
    T_IOSymbol out_symbol; //!< Output symbol corresponding to the edge
    float metric; //!< Metric of the edge
    CC_TreeNode<T_IOSymbol, T_Register, T_EdgeTag> *p_origin; //!< Pointer to the edge origin node
    CC_TreeNode<T_IOSymbol, T_Register, T_EdgeTag> *p_destination; //!< Pointer to the edge destination node
    T_EdgeTag edge_tag; //!< Optional and versatile class to tag the edge
};

} // namespace ccsoft


#endif // __CC_TREE_EDGE_H__

