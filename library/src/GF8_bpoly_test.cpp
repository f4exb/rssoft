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

	 Tests of Galois Field sub-library on GF(8)[X,Y]

*/

#include <iostream>
#include <utility>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include "GFq.h"
#include "GFq_Element.h"
#include "GFq_Polynomial.h"
#include "GF2_Element.h"
#include "GF2_Polynomial.h"
#include "GFq_BivariateMonomial.h"
#include "GFq_BivariatePolynomial.h"

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
	unsigned int k = 3;

	rssoft::gf::GFq_BivariateMonomial m0(rssoft::gf::GFq_Element(gf8,1),1,1); // X*Y
	rssoft::gf::GFq_BivariateMonomial m1(rssoft::gf::GFq_Element(gf8,1),2,0); // X^2
	rssoft::gf::GFq_BivariateMonomial m2(rssoft::gf::GFq_Element(gf8,1),0,0); // 1
	rssoft::gf::GFq_BivariateMonomial m3(rssoft::gf::GFq_Element(gf8,1),3,0); // a*X^3

	std::cout << m0 << " wdeg = " << m0.wdeg(1,k-1) << std::endl;
	std::cout << m1 << " wdeg = " << m1.wdeg(1,k-1) << std::endl;
	std::cout << m2 << " wdeg = " << m2.wdeg(1,k-1) << std::endl;
	std::cout << m3 << " wdeg = " << m3.wdeg(1,k-1) << std::endl;

	std::vector<rssoft::gf::GFq_BivariateMonomial> monos_P;
	monos_P.push_back(m0); // X*Y
	monos_P.push_back(m1); // X^2
	monos_P.push_back(m2); // 1
	monos_P.push_back(m3); // a*X^3

	rssoft::gf::GFq_BivariatePolynomial P(1,k-1);
	P.init(monos_P);

	std::cout << "P(X,Y) = " << P << std::endl;
}


