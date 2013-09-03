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

 Convolutional encoder class.

 This version uses a fixed array for internal registers

 */
#ifndef __CC_ENCODING_REGISTERS_FA_H__
#define __CC_ENCODING_REGISTERS_FA_H__

#include "CCSoft_Exception.h"

#include <iostream>
#include <array>
#include <algorithm>

namespace ccsoft
{

/**
 * \brief Convolutional encoding registers class.
 * \tparam T_Register type of the internal registers
 * \tparam N_k Size of an input symbol in bits (k parameter)
 */
template<typename T_Register, unsigned int N_k>
class CC_EncodingRegisters_FA
{
public:

    //=============================================================================================
    /**
     * Constructor.
     */
    CC_EncodingRegisters_FA()
    {
        clear();
    }

    //=============================================================================================
    /**
     * Destructor
     */
    ~CC_EncodingRegisters_FA()
    {}

    //=============================================================================================
    /**
     * Clear internal registers. Used before encoding a sequence
     */
    void clear()
    {
        std::fill(registers.begin(), registers.end(), 0);
        //registers.fill(0);
    }
    
    //=============================================================================================
    /**
     * Get a register given its index
     */
    T_Register& get_register(const unsigned int index)
    {
        return registers[index];
    }
    
    //=============================================================================================
    /**
     * Get registers reference
     */
    const std::array<T_Register, N_k>& get_registers() const
    {
        return registers;
    }

    //=============================================================================================
    /**
     * Set registers
     */
    void set_registers(const std::array<T_Register, N_k>& _registers)
    {
        registers = _registers;
    }


protected:
    std::array<T_Register, N_k> registers; //!< Memory registers as many as there are inputs
};

/**
 * \brief Convolutional encoding registers class. Specialized in single register.
 * \tparam T_Register type of the internal registers
 */
template<typename T_Register>
class CC_EncodingRegisters_FA<T_Register, 1>
{
public:

    //=============================================================================================
    /**
     * Constructor.
     */
    CC_EncodingRegisters_FA()
    {
        clear();
    }

    //=============================================================================================
    /**
     * Destructor
     */
    ~CC_EncodingRegisters_FA()
    {}

    //=============================================================================================
    /**
     * Clear internal registers. Used before encoding a sequence
     */
    void clear()
    {
        m_register = 0;
    }
    
    //=============================================================================================
    /**
     * Get a register given its index
     */
    T_Register& get_register(const unsigned int index)
    {
        if (index > 0)
        {
            throw CCSoft_Exception("Invalid subscript for single register");
        }
    
        return m_register;
    }
    
    //=============================================================================================
    /**
     * Get registers reference
     */
    const T_Register& get_registers() const
    {
        return m_register;
    }

    //=============================================================================================
    /**
     * Set registers
     */
    void set_registers(const T_Register& _register)
    {
        m_register = _register;
    }


protected:
    T_Register m_register; //!< Memory single register
};


} // namespace ccsoft


#endif // __CC_ENCODING_REGISTERS_FA_H__
