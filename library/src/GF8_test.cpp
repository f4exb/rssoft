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

	 Tests of Galois Field sub-library on GF(8) and GF(8)[X]

*/

#include <iostream>
#include <utility>
#include <stdlib.h>
#include <stdio.h>
#include "GFq.h"
#include "GFq_Element.h"
#include "GFq_Polynomial.h"
#include "GF2_Element.h"
#include "GF2_Polynomial.h"


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
template <typename T>
std::ostream& operator <<(std::ostream& os, const std::vector<T>& vec)
{
    os << "[";
    typename std::vector<T>::const_iterator it = vec.begin();

    for (;it != vec.end(); ++ it)
    {
    	os << (it == vec.begin() ? "" : ", ");
    	os << *it;
    }

    os << "]";

    return os;
}


// ================================================================================================
int main(int argc, char *argv[])
{

    std::cout << gf8;

    rssoft::gf::GFq_Element test(gf8), def_e(gf8), one(gf8,1);
    test = def_e^2;
    std::cout << test << std::endl;
    std::cout << one << std::endl;

    rssoft::gf::GFq_Element pe[4] = {
        rssoft::gf::GFq_Element(gf8,gf8.alpha(4)),
        rssoft::gf::GFq_Element(gf8,gf8.alpha(1)),
        rssoft::gf::GFq_Element(gf8,gf8.alpha(2)),
        rssoft::gf::GFq_Element(gf8,gf8.alpha(3)),
    };
    size_t pe_sz = sizeof(pe)/sizeof(rssoft::gf::GFq_Element);
    
    rssoft::gf::GFq_Element qe[3] = {
        rssoft::gf::GFq_Element(gf8,1),
        rssoft::gf::GFq_Element(gf8,0),
        rssoft::gf::GFq_Element(gf8,1),
    };
    size_t qe_sz = sizeof(qe)/sizeof(rssoft::gf::GFq_Element);
    
    rssoft::gf::GFq_Element qz[1] = {
        rssoft::gf::GFq_Element(gf8,1),
    };
    size_t qz_sz = sizeof(qz)/sizeof(rssoft::gf::GFq_Element);    
    
    std::cout << gf8.pwr() << std::endl;
    rssoft::gf::GFq_Polynomial P(gf8,pe_sz,pe);
    std::cout << "P(x) = " << P << " deg=" << P.deg() << std::endl;
    P.set_alpha_format(true);
    std::cout << "P(x) = " << P << std::endl;
    rssoft::gf::GFq_Polynomial Pc(P);
    std::cout << "Pc(x) = " << Pc << " deg=" << Pc.deg() << std::endl;
    Pc.set_alpha_format(false);
    std::cout << "Pc(x) = " << Pc << std::endl;
    rssoft::gf::GFq_Polynomial Q(gf8,qe_sz,qe);
    std::cout << "Q(x) = " << Q << " deg=" << Q.deg() << std::endl;
    Q.set_alpha_format(true);
    std::cout << "Q(x) = " << Q << std::endl;
   
    std::cout << "Q(a) = " << Q(gf8.alpha(1)) << std::endl;

    rssoft::gf::GFq_Polynomial S=P+Q;
    std::cout << "S(x) = " << S << std::endl;
    S.set_alpha_format(false);
    std::cout << "S(x) = " << S << std::endl;
   
    rssoft::gf::GFq_Polynomial M=P*Q;
    std::cout << "M(x) = " << M << std::endl;
    M.set_alpha_format(false);
    std::cout << "M(x) = " << M << std::endl;
   
    P.set_alpha_format(false);
    Q.set_alpha_format(false);
    const std::pair<rssoft::gf::GFq_Polynomial,rssoft::gf::GFq_Polynomial>& divres = div(P,Q);
    std::cout << "q(x) = " << divres.first << std::endl;
    std::cout << "r(x) = " << divres.second << std::endl;
    std::cout << "q(x)*Q(x)+r(x) = " << divres.first*Q+divres.second << std::endl;
   
    std::cout << "P(x)/Q(x) = " << P/Q << std::endl;
    std::cout << "P(x)%Q(x) = " << P%Q << std::endl;
    
    
    rssoft::gf::GFq_Polynomial Z(gf8,qz_sz,qz);
    std::cout << "Z(x) = " << Z << " deg=" << Z.deg() << std::endl;
    std::cout << "P(x)/Z(x) = " << P/Z << std::endl;
    
    rssoft::gf::GFq_Polynomial G=gcd(P,Q);
    std::cout << "gcd(P,Q)(x) = " << G << std::endl;

    rssoft::gf::GFq_Element ae[3] = {
        rssoft::gf::GFq_Element(gf8,0),
        rssoft::gf::GFq_Element(gf8,1),
        rssoft::gf::GFq_Element(gf8,1),
    };
    size_t ae_sz = sizeof(ae)/sizeof(rssoft::gf::GFq_Element);
    rssoft::gf::GFq_Polynomial A(gf8,ae_sz,ae);
    
    rssoft::gf::GFq_Element be[2] = {
        rssoft::gf::GFq_Element(gf8,1),
        rssoft::gf::GFq_Element(gf8,1),
    };
    size_t be_sz = sizeof(be)/sizeof(rssoft::gf::GFq_Element);
    rssoft::gf::GFq_Polynomial B(gf8,be_sz,be);

    std::cout << "gcd(A,B)(x) = " << gcd(A,B) << std::endl;
    
    std::cout << "P'(x) = " << P.derivative() << std::endl;
    
    std::cout << "P(x)<<1 = " << (P<<1) << std::endl;
    std::cout << "P(x)>>1 = " << (P>>1) << std::endl;
    std::cout << "P(x)<<2 = " << (P<<2) << std::endl;
    std::cout << "P(x)>>2 = " << (P>>2) << std::endl;
    std::cout << "P(x)<<3 = " << (P<<3) << std::endl;
    std::cout << "P(x)>>3 = " << (P>>3) << std::endl;
    
    const std::vector<rssoft::gf::GFq_Element>& rootsP = rootex(P);
    std::cout << "roots(P) = " << rootsP << std::endl;

    const std::vector<rssoft::gf::GFq_Element>& rootsQ = rootex(Q);
    std::cout << "roots(Q) = " << rootsQ << std::endl;

    const std::vector<rssoft::gf::GFq_Element>& rootsChienP = P.rootChien();
    std::cout << "(Chien's) roots(P) = " << rootsChienP << std::endl;

    const std::vector<rssoft::gf::GFq_Element>& rootsChienQ = Q.rootChien();
    std::cout << "(Chien's) roots(Q) = " << rootsChienQ << std::endl;

    rssoft::gf::GFq_Element ce[2] = {
        rssoft::gf::GFq_Element(gf8,6),
        rssoft::gf::GFq_Element(gf8,1),
    };
    size_t ce_sz = sizeof(ce)/sizeof(rssoft::gf::GFq_Element);
    rssoft::gf::GFq_Polynomial C(gf8,ce_sz,ce);

    const std::pair<rssoft::gf::GFq_Polynomial,rssoft::gf::GFq_Polynomial>& divpc = div(P,C);
    std::cout << "P(X) = (" << C << ")*(" << divpc.first << ") + (" << divpc.second << ")" << std::endl;

    const std::vector<rssoft::gf::GFq_Element>& rootsPC = rootex(divpc.first);
    std::cout << "roots(P/C) = " << rootsPC << std::endl;

    ce[0]=rssoft::gf::GFq_Element(gf8,1);
    rssoft::gf::GFq_Polynomial C1(gf8,ce_sz,ce);
    ce[0]=rssoft::gf::GFq_Element(gf8,2);
    rssoft::gf::GFq_Polynomial C2(gf8,ce_sz,ce);
    ce[0]=rssoft::gf::GFq_Element(gf8,3);
    rssoft::gf::GFq_Polynomial C3(gf8,ce_sz,ce);

    const std::vector<rssoft::gf::GFq_Element>& rootsMCi = rootex(C*C1*C2*C3);
    std::cout << "roots(C*C1*C2*C3) = " << rootsMCi << std::endl;

    const std::vector<rssoft::gf::GFq_Element>& rootsChienMCi = rootex(C*C1*C2*C3);
    std::cout << "(Chien's) roots(C*C1*C2*C3) = " << rootsChienMCi << std::endl;

    rssoft::gf::GFq_Polynomial D=P;
    const rssoft::gf::GFq_Element d_lead = D.make_monic();
    std::cout << "P.make_monic(X)) = " << D << " lead : " << d_lead << std::endl;
    rssoft::gf::GFq_Element d_lead2(D.field());
    D = get_monic(P, d_lead2);
    std::cout << "get_monic(P(X)) = " << D << " lead : " << d_lead2 << std::endl;
    
    std::cout << "P (D) square free decomposition: " << square_free_decomposition(D) << std::endl;
    rssoft::gf::GFq_Polynomial Cx = C*C1*C2*C3;
    Cx.make_monic();
    std::cout << "C*C1*C2*C3(X) (monic) = " << Cx << std::endl;
    std::cout << "C*C1*C2*C3 square free decomposition: " << square_free_decomposition(Cx) << std::endl;

    // use monic P
    /*
    const std::vector<rssoft::gf::GFq_Polynomial>& decomp_P = square_free_decomposition(D);
    std::vector<rssoft::gf::GFq_Polynomial>::const_iterator pd_it = decomp_P.begin();

    for (; pd_it != decomp_P.end(); ++pd_it)
    {
    	std::cout << "* " << *pd_it << std::endl;
    }
    */

    exit(EXIT_SUCCESS);
    return true;
}

