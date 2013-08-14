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

template<typename T_IOSymbol, typename T_Register>
class CC_TreeNode;

/**
 * \brief An edge of the code tree
 * \tparam T_IOSymbol Type of the input and output symbols
 * \tparam T_Register Type of the encoder internal registers
 */
template<typename T_IOSymbol, typename T_Register>
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
            CC_TreeNode<T_IOSymbol, T_Register> *_p_origin) :
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
    void set_p_destination(CC_TreeNode<T_IOSymbol, T_Register> *_p_destination)
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
    CC_TreeNode<T_IOSymbol, T_Register> *get_p_origin()
    {
        return p_origin;
    }

    /**
     * Destination pointer getter
     */
    CC_TreeNode<T_IOSymbol, T_Register> *get_p_destination()
    {
        return p_destination;
    }

protected:
    unsigned int id; //!< Edge's unique ID
    T_IOSymbol in_symbol; //!< Input symbol corresponding to the edge
    T_IOSymbol out_symbol; //!< Output symbol corresponding to the edge
    float metric; //!< Metric of the edge
    CC_TreeNode<T_IOSymbol, T_Register> *p_origin; //!< Pointer to the edge origin node
    CC_TreeNode<T_IOSymbol, T_Register> *p_destination; //!< Pointer to the edge destination node
};

} // namespace ccsoft


#endif // __CC_TREE_EDGE_H__

