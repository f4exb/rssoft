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

 Bivariate monomials with coefficient in GF(2^m) class

 */
#include "GFq_BivariateMonomial.h"

namespace rssoft
{
namespace gf
{

// ================================================================================================
GFq_BivariateMonomial::GFq_BivariateMonomial(const GFq_Element& coeff, unsigned int x_pow, unsigned int y_pow) :
		coefficient(coeff), exponents(x_pow, y_pow)
{}

// ================================================================================================
GFq_BivariateMonomial::~GFq_BivariateMonomial()
{}

// ================================================================================================
GFq_WeightedRevLex_BivariateMonomial::GFq_WeightedRevLex_BivariateMonomial(unsigned int w_x, unsigned int w_y) :
		weights(w_x, w_y)
{}

// ================================================================================================
GFq_WeightedRevLex_BivariateMonomial::~GFq_WeightedRevLex_BivariateMonomial()
{}

// ================================================================================================
bool GFq_WeightedRevLex_BivariateMonomial::operator()(const GFq_BivariateMonomial& m1, const GFq_BivariateMonomial& m2) const
{
	unsigned int w1 = m1.wdeg(weights.first, weights.second);
	unsigned int w2 = m2.wdeg(weights.first, weights.second);

	if (w1 == w2)
	{
		return m1.x_pow() > m2.x_pow();
	}
	else
	{
		return (w1 < w2);
	}
}

// ================================================================================================
std::ostream& operator <<(std::ostream& os, const GFq_BivariateMonomial& monomial)
{
	if (monomial.coefficient.is_zero())
	{
		if ((monomial.exponents.first == 0) && (monomial.exponents.second == 0))
		{
			os << "0";
		}
	}
	else if (monomial.coefficient == 1)
	{
		if ((monomial.exponents.first == 0) && (monomial.exponents.second == 0))
		{
			os << "1";
		}
	}
	else
	{
		os << monomial.coefficient << "*";
	}

	if (monomial.exponents.first > 0)
	{
		os << "X";

		if (monomial.exponents.first > 1)
		{
			os << "^" << monomial.exponents.first;
		}
	}

	if (monomial.exponents.second > 0)
	{
		os << "*Y";

		if (monomial.exponents.second > 1)
		{
			os << "^" << monomial.exponents.second;
		}
	}

	return os;
}

} // namespace gf
} // namsepace rssoft
