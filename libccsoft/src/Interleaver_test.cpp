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

     Interleaver test

*/

#include "CC_StackDecoding.h"
#include "CCSoft_Exception.h"

// ================================================================================================
// template to print a vector of printable elements
template<typename TDisplay, typename TElement, typename TStream> void print_vector(const std::vector<TElement>& v, TStream& os)
{
    os << "[";

    typename std::vector<TElement>::const_iterator it = v.begin();
    const typename std::vector<TElement>::const_iterator v_end = v.end();

    for (; it != v_end; ++it)
    {

        os <<  (TDisplay)(*it);

        if (it != v.begin()+v.size()-1)
        {
            os << ", ";
        }

    }

    os << "]";
}

int main(int argc, char *argv[])
{
    try
    {
        std::vector<unsigned int> k_constraints(1,32);
        std::vector<unsigned int> g1;
        g1.push_back(4073739089);
        g1.push_back(3831577671);
        std::vector<std::vector<unsigned int> > generator_polys(1,g1);

        ccsoft::CC_StackDecoding<unsigned int, unsigned int> cc_decoding(k_constraints, generator_polys);

        std::vector<unsigned int> test_symbols;

        for (unsigned int i=0; i<10; i++)
        {
            test_symbols.push_back(i);
        }

        print_vector<unsigned int, unsigned int, std::ostream>(test_symbols, std::cout);
        std::cout << std::endl;

        std::cout << "interleave:" << std::endl;
        cc_decoding.interleave(test_symbols);
        print_vector<unsigned int, unsigned int, std::ostream>(test_symbols, std::cout);
        std::cout << std::endl;

        std::cout << "de-interleave:" << std::endl;
        cc_decoding.interleave(test_symbols, false);
        print_vector<unsigned int, unsigned int, std::ostream>(test_symbols, std::cout);
        std::cout << std::endl;
    }
    catch (ccsoft::CCSoft_Exception& e)
    {
        std::cout << "CCSoft exception caught: " << e.what() << std::endl;
    }
}

