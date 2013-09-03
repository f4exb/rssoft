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

 Convolutional encoder class

 */
#ifndef __CC_ENCODING_H__
#define __CC_ENCODING_H__

#include "CC_Encoding_base.h"
#include "CCSoft_Exception.h"
#include <vector>
#include <iostream>

namespace ccsoft
{

/**
 * \brief Convolutional encoding class. Supports any k,n with k<n.
 * The input bits of a symbol are clocked simultaneously into the right hand side, or least significant position of the internal
 * registers. Therefore the given polynomial representation of generators should follow the same convention.
 * \tparam T_Register type of the internal registers
 * \tparam T_IOSymbol type used to pass input and output symbols
 */
template<typename T_Register, typename T_IOSymbol>
class CC_Encoding : public CC_Encoding_base<T_Register, T_IOSymbol>
{
public:

    //=============================================================================================
    /**
     * Constructor.
     * \param _constraints Vector of register lengths (constraint length + 1). The number of elements determines k.
     * \param _genpoly_representations Generator polynomial numeric representations. There are as many elements as there
     * are input bits (k). Each element is itself a vector with one polynomial value per output bit. The smallest size of
     * these vectors is retained as the number of output bits n. The input bits of a symbol are clocked simultaneously into
     * the right hand side, or least significant position of the internal registers. Therefore the given polynomial representation
     * of generators should follow the same convention.
     */
    CC_Encoding(const std::vector<unsigned int>& _constraints, const std::vector<std::vector<T_Register> >& _genpoly_representations) :
        CC_Encoding_base<T_Register, T_IOSymbol>(_constraints, _genpoly_representations),
        registers(_constraints.size(),0)
    {
    }

    //=============================================================================================
    /**
     * Destructor
     */
    virtual ~CC_Encoding()
    {}

    //=============================================================================================
    /**
     * Clear internal registers. Used before encoding a sequence
     */
    void clear()
    {
        registers.assign(registers.size(), 0);
    }

    //=============================================================================================
    /**
     * Get a R/W reference to a regiser
     * \param index Index of the register
     */
    virtual T_Register& get_register(unsigned int index)
    {
        return registers[index];
    }
    
    /**
     * Get registers reference
     */
    const std::vector<T_Register>& get_registers() const
    {
        return registers;
    }

    /**
     * Set registers
     */
    void set_registers(const std::vector<T_Register>& _registers)
    {
        registers = _registers;
    }


protected:
    std::vector<T_Register> registers; //!< Memory registers as many as there are inputs
};

} // namespace ccsoft


#endif // __CC_ENCODING_H__
