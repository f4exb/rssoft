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
#ifndef __CC_ENCODER_H__
#define __CC_ENCODER_H__

#include "CCSoft_Exception.h"
#include <vector>
#include <iostream>

namespace ccsoft
{

template<typename T_Register>
bool xorbits(const T_Register& reg)
{
    T_Register w_reg = reg;
    unsigned int nb_ones = 0;

    while(w_reg != 0)
    {
        nb_ones += w_reg % 2;
        w_reg /= 2;
    }

    return (nb_ones % 2) == 1;
}

template<typename T_Register>
void print_register(const T_Register& reg, std::ostream& os)
{
    os << std::hex << reg;
}

template<>
void print_register<unsigned char>(const unsigned char& reg, std::ostream& os)
{
    os << std::hex << (unsigned int) reg;
}

/**
 * \brief Convolutional encoding class. Supports any k,n with k<n.
 * The input bits of a symbol are clocked simultaneously into the right hand side, or least significant position of the internal
 * registers. Therefore the given polynomial representation of generators should follow the same convention.
 * \tparam T_Register type of the internal registers
 * \tparam T_IOSymbol type used to pass input and output symbols
 */
template<typename T_Register, typename T_IOSymbol>
class CC_Encoding
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
        k(_constraints.size()),
        constraints(_constraints),
        genpoly_representations(_genpoly_representations),
        registers(k,0)
    {
        if (k < 1)
        {
            throw CCSoft_Exception("There must be at least one constraint size");
        }

        if (k > sizeof(T_IOSymbol)*8)
        {
            throw CCSoft_Exception("Number of input bits not supported by I/O symbol type");
        }

        if (genpoly_representations.size() != k)
        {
            throw CCSoft_Exception("Generator polynomial representations size error");
        }

        unsigned int min_nb_outputs = genpoly_representations[0].size();

        for (unsigned int ci=0; ci < constraints.size(); ci++)
        {
            if (constraints[ci] > sizeof(T_Register)*8)
            {
                throw CCSoft_Exception("One constraint size is too large for the size of the registers");
            }

            if (genpoly_representations[ci].size() < min_nb_outputs)
            {
                min_nb_outputs = genpoly_representations[ci].size();
            }
        }

        n = min_nb_outputs;

        if (n <= k)
        {
            throw CCSoft_Exception("The number of outputs must be larger than the number of inputs");
        }

        if (n > sizeof(T_IOSymbol)*8)
        {
            throw CCSoft_Exception("Number of output bits not supported by I/O symbol type");
        }
    }

    //=============================================================================================
    /**
     * Destructor
     */
    ~CC_Encoding()
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
     * Encode a new symbol of k bits into a symbol of n bits
     * \param in_symbol Input symbol
     * \param out_symbol Output symbol
     */
    bool encode(const T_IOSymbol& in_symbol, T_IOSymbol& out_symbol)
    {
        T_IOSymbol w_in = in_symbol;

        // load registers with new symbol bits
        for (unsigned int ki=0; ki<k; ki++)
        {
            registers[ki] <<= 1; // make room for bit
            registers[ki] += w_in & 1; // insert bit
            w_in >>= 1; // move to next input bit
        }

        out_symbol = 0;
        T_IOSymbol symbol_bit;

        for (unsigned int ni=0; ni<n; ni++)
        {
            symbol_bit = 0;

            for (unsigned int ki=0; ki<k; ki++)
            {
                symbol_bit ^= (xorbits(registers[ki]&genpoly_representations[ki][ni]) ? 1 : 0);
            }

            out_symbol += symbol_bit << ni;
        }

        return true;
    }

    //=============================================================================================
    /**
     * Prints encoding characteristics to an output stream
     */
    void print(std::ostream& os)
    {
        std::cout << "k=" << k << ", n=" << n << std::endl;

        for (unsigned int ci=0; ci<k; ci++)
        {
            os << ci << " (" << constraints[ci] << ") : ";

            for (unsigned int gi=0; gi<n; gi++)
            {
                print_register(genpoly_representations[ci][gi], os);
                os << " ";
            }

            os << std::endl;
        }
    }

protected:
    unsigned int k;
    unsigned int n;
    std::vector<unsigned int> constraints; //!< As many constraints as there are inputs
    std::vector<std::vector<T_Register> > genpoly_representations; //!< As many generator polynomials vectors (the size of the number of outputs) as there are inputs
    std::vector<T_Register> registers; // Memory registers as many as there are inputs
};

} // namespace ccsoft


#endif // __CC_ENCODER_H__
