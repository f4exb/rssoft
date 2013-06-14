/*
     Copyright 2013 Edouard Griffiths <f4exb at free dot fr>

     This file is part of RSSoft. A Reed-Solomon Soft Decoding library

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

	 Original from Arash Partow (see original copytight notice).
	 Modified and included in RSSoft in gf sub-namespace.

	 Unit test of decoder based on the Soft decision decoding complete example
	 in the Reed-Solomon codes Part 2 document

*/

#include "GF2_Element.h"
#include "GF2_Polynomial.h"
#include "GFq.h"
#include "ReliabilityMatrix.h"
#include "MultiplicityMatrix.h"
#include "GSKV_Interpolation.h"
#include <iostream>
#include <iomanip>

float pwr_S0[8] = {2.163577, 0.003943, 0.064378, 0.000117, 0.021512, 0.000038, 0.000640, 0.000002};
float pwr_S1[8] = {0.459689, 0.012363, 0.011172, 0.000300, 1.580876, 0.042520, 0.038420, 0.001032};
float pwr_S2[8] = {0.009034, 0.000000, 0.000245, 0.000000, 1.603912, 0.000010, 0.043565, 0.000000};
float pwr_S3[8] = {0.736172, 0.838307, 0.005258, 0.005987, 0.005077, 0.005782, 0.000037, 0.000042};
float pwr_S4[8] = {0.001144, 0.912521, 0.000128, 0.102537, 0.000000, 0.000312, 0.000000, 0.000036};
float pwr_S5[8] = {0.000708, 0.036403, 0.026054, 1.339624, 0.000000, 0.000004, 0.000003, 0.000129};
float pwr_S6[8] = {1.507900, 0.000456, 0.045338, 0.000013, 0.607732, 0.000183, 0.018272, 0.000007};

/*
   P(X) = X^3+X+1
   p(x) = 1x^3+0x^2+1x^1+1x^0
             1    0    1    1 <-MSB
*/
rssoft::gf::GF2_Element ppe[4] = {1,1,0,1};
rssoft::gf::GF2_Polynomial ppoly(4,ppe);
/*
  A Galois Field of type GF(2^3)
*/
rssoft::gf::GFq gf8(3,ppoly);

// ================================================================================================
int main(int argc, char *argv[])
{
	rssoft::ReliabilityMatrix mat_Pi(3,7);

	std::cout << "Power matrix:" << std::endl;
	mat_Pi.enter_symbol_data(pwr_S0);
	mat_Pi.enter_symbol_data(pwr_S1);
	mat_Pi.enter_symbol_data(pwr_S2);
	mat_Pi.enter_symbol_data(pwr_S3);
	mat_Pi.enter_symbol_data(pwr_S4);
	mat_Pi.enter_symbol_data(pwr_S5);
	mat_Pi.enter_symbol_data(pwr_S6);

    std::cout << mat_Pi;

	mat_Pi.normalize();

	std::cout << std::endl;
	std::cout << "Reliability matrix:" << std::endl;
    std::cout << mat_Pi;

    rssoft::MultiplicityMatrix mat_M_f(mat_Pi, 3.0f);
    
	std::cout << std::endl;
	std::cout << "Multiplicity matrix (short construction):" << std::endl;
    std::cout << mat_M_f;
    std::cout << "cost is " << mat_M_f.cost() << std::endl;
    
    rssoft::MultiplicityMatrix mat_M(mat_Pi, 12u);

	std::cout << std::endl;
	std::cout << "Multiplicity matrix (long construction):" << std::endl;
    std::cout << mat_M;
    std::cout << "cost is " << mat_M.cost() << std::endl;
    
    std::cout << std::endl;
    rssoft::MultiplicityMatrix::traversing_iterator m_it(mat_M.begin());
    for (; m_it != mat_M.end(); ++m_it)
    {
        std::cout << "M(" << m_it.iX() << "," << m_it.iY() << ") = " << m_it.multiplicity() << std::endl;
    }

    rssoft::GSKV_Interpolation gskv(gf8, 5);
    gskv.run(mat_M);
}
