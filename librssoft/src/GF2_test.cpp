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

	 Tests of Galois Field sub-library on GF(2)[X]

*/

#include <iostream>
#include <iomanip>
#include <utility>
#include <stdlib.h>
#include <stdio.h>
#include "GF2_Polynomial.h"

void show_irred_prim(rssoft::gf::GF2_Polynomial& a);

int main(int argc, char *argv[])
{

    rssoft::gf::GF2_Element pe[4] = {
        rssoft::gf::GF2_Element(1),
        rssoft::gf::GF2_Element(0),
        rssoft::gf::GF2_Element(1),
        rssoft::gf::GF2_Element(1),
    };
    size_t pe_sz = sizeof(pe)/sizeof(rssoft::gf::GF2_Element);
    
    rssoft::gf::GF2_Polynomial P(pe_sz,pe); // P(X) = 1 + x^2 + x^3
    std::cout << "P(x) = " << P << " deg=" << P.deg() << std::endl;

    rssoft::gf::GF2_Element qe[3] = {
        rssoft::gf::GF2_Element(1),
        rssoft::gf::GF2_Element(0),
        rssoft::gf::GF2_Element(1),
    };
    size_t qe_sz = sizeof(qe)/sizeof(rssoft::gf::GF2_Element); // Q(X) = 1 + x^2

    rssoft::gf::GF2_Polynomial Q(qe_sz,qe);
    std::cout << "Q(x) = " << Q << " deg=" << Q.deg() << std::endl;

    rssoft::gf::GF2_Element se[4] = {
        rssoft::gf::GF2_Element(1),
        rssoft::gf::GF2_Element(1),
        rssoft::gf::GF2_Element(0),
        rssoft::gf::GF2_Element(1),
    };
    size_t se_sz = sizeof(se)/sizeof(rssoft::gf::GF2_Element); // Q(X) = 1 + x + x^3

    rssoft::gf::GF2_Polynomial S(se_sz,se);
    std::cout << "S(x) = " << S << " deg=" << S.deg() << std::endl;

    rssoft::gf::GF2_Element te[4] = {
        rssoft::gf::GF2_Element(1),
        rssoft::gf::GF2_Element(0),
        rssoft::gf::GF2_Element(0),
        rssoft::gf::GF2_Element(1),
    };
    size_t te_sz = sizeof(te)/sizeof(rssoft::gf::GF2_Element); // T(X) = 1 + x^3

    rssoft::gf::GF2_Polynomial T(te_sz,te);
    std::cout << "T(x) = " << T << " deg=" << T.deg() << std::endl;

    rssoft::gf::GF2_Polynomial M = P*Q;
    std::cout << "M(x) = " << M << std::endl;

    const std::pair<rssoft::gf::GF2_Polynomial,rssoft::gf::GF2_Polynomial>& divres = div(M,Q);
    std::cout << "q(x) = " << divres.first << std::endl;
    std::cout << "r(x) = " << divres.second << std::endl;
	
	std::cout << "gcd(Q*P,Q*T)(x) = " << gcd(Q*P,Q*T) << std::endl;
	std::cout << "gcd(Q*P,Q*P)(x) = " << gcd(Q*P,Q*P) << std::endl;
	
    rssoft::gf::GF2_Element ue[6] = {
        rssoft::gf::GF2_Element(1),
        rssoft::gf::GF2_Element(1),
        rssoft::gf::GF2_Element(1),
        rssoft::gf::GF2_Element(0),
        rssoft::gf::GF2_Element(1),
        rssoft::gf::GF2_Element(1),
    };
    size_t ue_sz = sizeof(ue)/sizeof(rssoft::gf::GF2_Element); // U(X) = 1 + x + x^2 + x^4 + x^5

    rssoft::gf::GF2_Polynomial U(ue_sz,ue);
    std::cout << "U(x) = " << U << " deg=" << U.deg() << std::endl;
	
    rssoft::gf::GF2_Element ve[3] = {
        rssoft::gf::GF2_Element(1),
        rssoft::gf::GF2_Element(1),
        rssoft::gf::GF2_Element(1),
    };
    size_t ve_sz = sizeof(ve)/sizeof(rssoft::gf::GF2_Element); // V(X) = 1 + x + x^2

    rssoft::gf::GF2_Polynomial V(ve_sz,ve);
    std::cout << "V(x) = " << V << " deg=" << V.deg() << std::endl;
	
	rssoft::gf::GF2_Element Xe[2] = {
		rssoft::gf::GF2_Element(0),
		rssoft::gf::GF2_Element(1),
	};
	size_t Xe_sz = sizeof(Xe)/sizeof(rssoft::gf::GF2_Element);
	rssoft::gf::GF2_Polynomial X(Xe_sz,Xe); // X(x)=x
	
	rssoft::gf::GF2_Polynomial A=U-X;
	
    rssoft::gf::GF2_Element be[5] = {
        rssoft::gf::GF2_Element(1),
        rssoft::gf::GF2_Element(1),
        rssoft::gf::GF2_Element(1),
        rssoft::gf::GF2_Element(1),
        rssoft::gf::GF2_Element(1),
    };
    size_t be_sz = sizeof(be)/sizeof(rssoft::gf::GF2_Element); 

    rssoft::gf::GF2_Polynomial B(be_sz,be);
	
	show_irred_prim(P);
	show_irred_prim(Q);
	show_irred_prim(S);
	show_irred_prim(T);
	show_irred_prim(U);
	show_irred_prim(A);
	show_irred_prim(B);
	show_irred_prim(V);
	
	/*
	std::cout << "is P(X) irreducible? " << (irreducible(P) ? "yes" : "no") << std::endl;
	std::cout << "is Q(X) irreducible? " << (irreducible(Q) ? "yes" : "no") << std::endl;
	std::cout << "is S(X) irreducible? " << (irreducible(S) ? "yes" : "no") << std::endl;
	std::cout << "is T(X) irreducible? " << (irreducible(T) ? "yes" : "no") << std::endl;
	std::cout << "is U(X) irreducible? " << (irreducible(U) ? "yes" : "no") << std::endl;
	std::cout << "is U(X)-X irreducible? " << (irreducible(U-X) ? "yes" : "no") << std::endl;
	std::cout << "is V(X) irreducible? " << (irreducible(V) ? "yes" : "no") << std::endl;
	*/

    exit(EXIT_SUCCESS);
    return true;

}

void show_irred_prim(rssoft::gf::GF2_Polynomial& P)
{
	std::cout << P << " -> (" << (irreducible(P) ? "Y" : "N") << "," << (primitive(P,P.deg()) ? "Y" : "N") << ")" << std::endl;
}
