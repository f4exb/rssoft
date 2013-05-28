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

 Bivariate polynomials with coefficients in GF(2^m) class

 */

#include "GFq_BivariatePolynomial.h"

namespace rssoft
{
namespace gf
{

// ================================================================================================
GFq_BivariatePolynomial::GFq_BivariatePolynomial(unsigned int w_x, unsigned int w_y) :
		weights(w_x,w_y),
		monomials(GFq_WeightedRevLex_BivariateMonomial(w_x,w_y))
{}

// ================================================================================================
GFq_BivariatePolynomial::~GFq_BivariatePolynomial()
{}

// ================================================================================================
void GFq_BivariatePolynomial::init(std::vector<GFq_BivariateMonomial>& _monomials)
{
	monomials.clear();
	monomials.insert(_monomials.begin(), _monomials.end());
}

// ================================================================================================
std::ostream& operator <<(std::ostream& os, const GFq_BivariatePolynomial& polynomial)
{
	if (polynomial.monomials.size() == 0)
	{
		os << "<invalid>";
	}
	else
	{
		std::set<GFq_BivariateMonomial>::const_iterator it = polynomial.monomials.begin();

		for (; it != polynomial.monomials.end(); ++it)
		{
			if (it != polynomial.monomials.begin())
			{
				os << " + ";
			}

			os << *it;
		}
	}

	return os;
}

} // namespace gf
} // namespace rssoft
