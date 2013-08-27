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

     Encoder tests

*/

#include "CC_StackDecoding.h"
#include "CCSoft_Exception.h"
#include "CC_ReliabilityMatrix.h"
#include <iostream>
#include <sstream>
#include <fstream>

int main(int argc, char *argv[])
{
    try
    {
        // Yunghsiang S. Han and Po-Ning Chen Sequential Decoding of Convolutional Codes
        // http://web.ntpu.edu.tw/~yshan/book_chapter.pdf
        // fig 1 example (2,1,2) code

        std::vector<unsigned int> hanchen1_ks(1,3);
        std::vector<unsigned char> hanchen1_g;
        hanchen1_g.push_back(7);
        hanchen1_g.push_back(5);
        std::vector<std::vector<unsigned char> > hanchen1_gs(1, hanchen1_g);

        ccsoft::CC_StackDecoding<unsigned char, unsigned char> hanchen1_cc_decoder(hanchen1_ks, hanchen1_gs);

        std::cout << "Han & Chen example 1:" << std::endl;
        hanchen1_cc_decoder.get_encoding().print(std::cout);

        unsigned char hanchen1_in_symbols[7] = {1,1,1,0,1,0,0};
        float hanchen1_soft_array[4] = {0.1,0.1,0.1,0.1};
        unsigned char out_symbol;
        std::ostringstream oos;
        ccsoft::CC_ReliabilityMatrix relmat_hanchen1_1(2, 7);

        for (unsigned int i=0; i<7; i++)
        {
            hanchen1_cc_decoder.get_encoding().encode(hanchen1_in_symbols[i], out_symbol);
            hanchen1_soft_array[out_symbol] = 0.7;
            relmat_hanchen1_1.enter_symbol_data(hanchen1_soft_array);
            hanchen1_soft_array[out_symbol] = 0.1;
            std::cout << (unsigned int) hanchen1_in_symbols[i] << " ";
            oos << (unsigned int) out_symbol << " ";
        }

        std::cout << std::endl;
        std::cout << oos.str() << std::endl;
        std::vector<unsigned char> result;
        hanchen1_cc_decoder.decode(relmat_hanchen1_1, result); // easy

        std::vector<unsigned char>::const_iterator r_it = result.begin();
        for (; r_it != result.end(); ++r_it)
        {
            std::cout << (unsigned int) *r_it << " ";
        }

        std::cout << std::endl << std::endl;

        std::ofstream hanchen1_dot_file;
        hanchen1_dot_file.open("hanchen1.dot");
        hanchen1_cc_decoder.print_dot(hanchen1_dot_file);
        hanchen1_dot_file.close();

        // fig 2 example (3,2,2) systematic code

        std::vector<unsigned int> hanchen2_ks(2,3);

        std::vector<unsigned char> hanchen2_g0;
        hanchen2_g0.push_back(1);
        hanchen2_g0.push_back(0);
        hanchen2_g0.push_back(2);

        std::vector<unsigned char> hanchen2_g1;
        hanchen2_g1.push_back(0);
        hanchen2_g1.push_back(1);
        hanchen2_g1.push_back(6);

        std::vector<std::vector<unsigned char> > hanchen2_gs;
        hanchen2_gs.push_back(hanchen2_g0);
        hanchen2_gs.push_back(hanchen2_g1);

        ccsoft::CC_StackDecoding<unsigned char, unsigned char> hanchen2_cc_decoder(hanchen2_ks, hanchen2_gs);

        std::cout << "Han & Chen example 2:" << std::endl;
        hanchen2_cc_decoder.get_encoding().print(std::cout);

        unsigned char hanchen2_in_symbols[4] = {3,2,0,0};
        float hanchen2_soft_array[8] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
        oos.str("");
        ccsoft::CC_ReliabilityMatrix relmat_hanchen2_1(3, 4);

        for (unsigned int i=0; i<4; i++)
        {
            hanchen2_cc_decoder.get_encoding().encode(hanchen2_in_symbols[i], out_symbol);
            hanchen2_soft_array[out_symbol] = 0.3;
            relmat_hanchen2_1.enter_symbol_data(hanchen2_soft_array);
            hanchen2_soft_array[out_symbol] = 0.1;
            std::cout << (unsigned int) hanchen2_in_symbols[i] << " ";
            oos << (unsigned int) out_symbol << " ";
        }

        std::cout << std::endl;
        std::cout << oos.str() << std::endl;

        result.clear();
        hanchen2_cc_decoder.decode(relmat_hanchen2_1, result); // easy

        r_it = result.begin();
        for (; r_it != result.end(); ++r_it)
        {
            std::cout << (unsigned int) *r_it << " ";
        }

        std::cout << std::endl << std::endl;

        std::ofstream hanchen2_dot_file;
        hanchen2_dot_file.open("hanchen2.dot");
        hanchen2_cc_decoder.print_dot(hanchen2_dot_file);
        hanchen2_dot_file.close();
    }
    catch (ccsoft::CCSoft_Exception& e)
    {
        std::cout << "CCSoft exception caught: " << e.what() << std::endl;
    }
}
