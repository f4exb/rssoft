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

     Full coding/decoding test with direct AWGN

*/

#include "CC_Encoding.h"
#include "CCSoft_Exception.h"
#include <iostream>

int main(int argc, char *argv[])
{
    try
    {
        // WSJT Layland-Lusbaugh constraint length 32 k=1 n=2

        std::vector<unsigned int> jt_ks(1,32);
        std::vector<unsigned int> jt_g;
        jt_g.push_back(0xf2d05351); // Layland-Lusbaugh G0
        jt_g.push_back(0xe4613c47); // Layland-Lusbaugh G1
        std::vector<std::vector<unsigned int> > jt_gs(1, jt_g);

        ccsoft::CC_Encoding<unsigned int, unsigned char> jt_cc_encoder(jt_ks, jt_gs);

        std::cout << "JT CC encoder:" << std::endl;
        jt_cc_encoder.print(std::cout);
        std::cout << std::endl;

        // Example given for Mathworks poly2treillis function

        std::vector<unsigned int> p2t_ks;
        p2t_ks.push_back(5);
        p2t_ks.push_back(4);

        std::vector<unsigned char> p2t_g0;
        p2t_g0.push_back(23);
        p2t_g0.push_back(35);
        p2t_g0.push_back(0);

        std::vector<unsigned char> p2t_g1;
        p2t_g1.push_back(0);
        p2t_g1.push_back(5);
        p2t_g1.push_back(13);

        std::vector<std::vector<unsigned char> > p2t_gs;
        p2t_gs.push_back(p2t_g0);
        p2t_gs.push_back(p2t_g1);

        ccsoft::CC_Encoding<unsigned char, unsigned char> p2t_cc_encoder(p2t_ks, p2t_gs);

        std::cout << "Mathworks poly2treillis example:" << std::endl;
        p2t_cc_encoder.print(std::cout);
        std::cout << std::endl;

        // Yunghsiang S. Han and Po-Ning Chen Sequential Decoding of Convolutional Codes
        // http://web.ntpu.edu.tw/~yshan/book_chapter.pdf
        // fig 1 example (2,1,2) code

        std::vector<unsigned int> hanchen1_ks(1,3);
        std::vector<unsigned char> hanchen1_g;
        hanchen1_g.push_back(7);
        hanchen1_g.push_back(5);
        std::vector<std::vector<unsigned char> > hanchen1_gs(1, hanchen1_g);

        ccsoft::CC_Encoding<unsigned char, unsigned char> hanchen1_cc_encoder(hanchen1_ks, hanchen1_gs);

        std::cout << "Han & Chen example 1:" << std::endl;
        hanchen1_cc_encoder.print(std::cout);

        unsigned char hanchen1_in_symbols[7] = {1,1,1,0,1,0,0};
        unsigned char out_symbol;

        for (unsigned int i=0; i<7; i++)
        {
            hanchen1_cc_encoder.encode(hanchen1_in_symbols[i], out_symbol);
            std::cout << (unsigned int) out_symbol << " ";
        }

        std::cout << std::endl << std::endl;

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

        ccsoft::CC_Encoding<unsigned char, unsigned char> hanchen2_cc_encoder(hanchen2_ks, hanchen2_gs);

        std::cout << "Han & Chen example 2:" << std::endl;
        hanchen2_cc_encoder.print(std::cout);

        unsigned char hanchen2_in_symbols[4] = {3,2,0,0};

        for (unsigned int i=0; i<4; i++)
        {
            hanchen2_cc_encoder.encode(hanchen2_in_symbols[i], out_symbol);
            std::cout << (unsigned int) out_symbol << " ";
        }

        std::cout << std::endl;
    }
    catch (ccsoft::CCSoft_Exception& e)
    {
        std::cout << "CCSoft exception caught: " << e.what() << std::endl;
    }
}
