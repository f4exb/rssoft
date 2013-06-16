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
GFq_BivariateMonomialExponents::GFq_BivariateMonomialExponents(unsigned int eX, unsigned int eY) :
	std::pair<unsigned int, unsigned int>(eX, eY)
{}

// ================================================================================================
GFq_BivariateMonomialExponents::GFq_BivariateMonomialExponents(const std::pair<unsigned int, unsigned int>& _exponents) :
	std::pair<unsigned int, unsigned int>(_exponents)
{}

// ================================================================================================
GFq_BivariateMonomial::GFq_BivariateMonomial(const GFq_Element& coeff, unsigned int eX, unsigned int eY) :
	std::pair<GFq_BivariateMonomialExponents, GFq_Element>(GFq_BivariateMonomialExponents(eX, eY), coeff)
{}

// ================================================================================================
GFq_BivariateMonomial::GFq_BivariateMonomial(const GFq_Element& coeff, const GFq_BivariateMonomialExponents& exponents) :
	std::pair<GFq_BivariateMonomialExponents, GFq_Element>(exponents, coeff)
{}

// ================================================================================================
GFq_BivariateMonomial::GFq_BivariateMonomial(const GFq_BivariateMonomialKeyValueRepresentation& monomial_rep) :
	std::pair<GFq_BivariateMonomialExponents, GFq_Element>(monomial_rep)
{}

// ================================================================================================
GFq_BivariateMonomial& GFq_BivariateMonomial::operator =(const GFq_BivariateMonomial& monomial)
{
	*this = monomial;
	return *this;
}

// ================================================================================================
GFq_BivariateMonomial& GFq_BivariateMonomial::operator +=(const GFq_BivariateMonomial& monomial)
{
	if (first != monomial.first)
	{
		throw GF_Exception("Cannot add monomials of different exponents");
	}
	else
	{
		second += monomial.second;
	}

	return *this;
}

// ================================================================================================
GFq_BivariateMonomial& GFq_BivariateMonomial::operator+=(const GFq_Element& gfe)
{
	second += gfe;
	return *this;
}

// ================================================================================================
GFq_BivariateMonomial& GFq_BivariateMonomial::operator -=(const GFq_BivariateMonomial& monomial)
{
	*this += monomial;
	return *this;
}

// ================================================================================================
GFq_BivariateMonomial& GFq_BivariateMonomial::operator-=(const GFq_Element& gfe)
{
	*this += gfe;
	return *this;
}

// ================================================================================================
GFq_BivariateMonomial& GFq_BivariateMonomial::operator *=(const GFq_BivariateMonomial& monomial)
{
	second *= monomial.second;
	first.first += monomial.first.first;
	first.second += monomial.first.second;
	return *this;
}

// ================================================================================================
GFq_BivariateMonomial& GFq_BivariateMonomial::operator*=(const GFq_Element& gfe)
{
	second *= gfe;
	return *this;
}

// ================================================================================================
GFq_BivariateMonomial& GFq_BivariateMonomial::operator /=(const GFq_BivariateMonomial& monomial)
{
	if (monomial.second == 0)
	{
		throw GF_Exception("Zero divide monomial");
	}

	second /= monomial.second;

	if (first.first - monomial.first.first < 0)
	{
		throw GF_Exception("Cannot divide by a monomial with a higher degree in X");
	}
	else if (first.second - monomial.first.second < 0)
	{
		throw GF_Exception("Cannot divide by a monomial with a higher degree in Y");
	}
	else
	{
		first.first -= monomial.first.first;
		first.second -= monomial.first.second;
	}

	return *this;
}

// ================================================================================================
GFq_BivariateMonomial& GFq_BivariateMonomial::operator/=(const GFq_Element& gfe)
{
	if (gfe == 0)
	{
		throw GF_Exception("Zero divide monomial");
	}

	second /= gfe;
	return *this;
}

// ================================================================================================
GFq_BivariateMonomial operator +(const GFq_BivariateMonomial& a, const GFq_BivariateMonomial& b)
{
	GFq_BivariateMonomial result(a);
	result += b;
	return result;
}

// ================================================================================================
GFq_BivariateMonomial operator +(const GFq_BivariateMonomial& a, const GFq_Element& b)
{
	GFq_BivariateMonomial result(a);
	result += b;
	return result;
}

// ================================================================================================
GFq_BivariateMonomial operator +(const GFq_Element& a, const GFq_BivariateMonomial& b)
{
	GFq_BivariateMonomial result(b);
	result += a;
	return result;
}

// ================================================================================================
GFq_BivariateMonomial operator -(const GFq_BivariateMonomial& a, const GFq_BivariateMonomial& b)
{
	GFq_BivariateMonomial result(a);
	result -= b;
	return result;
}

// ================================================================================================
GFq_BivariateMonomial operator -(const GFq_BivariateMonomial& a, const GFq_Element& b)
{
	GFq_BivariateMonomial result(a);
	result -= b;
	return result;
}

// ================================================================================================
GFq_BivariateMonomial operator -(const GFq_Element& a, const GFq_BivariateMonomial& b)
{
	GFq_BivariateMonomial result(b);
	result -= a;
	return result;
}

// ================================================================================================
GFq_BivariateMonomial operator *(const GFq_BivariateMonomial& a, const GFq_BivariateMonomial& b)
{
	GFq_BivariateMonomial result(a);
	result *= b;
	return result;
}

// ================================================================================================
GFq_BivariateMonomial operator *(const GFq_BivariateMonomial& a, const GFq_Element& b)
{
	GFq_BivariateMonomial result(a);
	result *= b;
	return result;
}

// ================================================================================================
GFq_BivariateMonomial operator *(const GFq_Element& a, const GFq_BivariateMonomial& b)
{
	GFq_BivariateMonomial result(b);
	result *= a;
	return result;
}

// ================================================================================================
GFq_BivariateMonomial operator /(const GFq_BivariateMonomial& a, const GFq_BivariateMonomial& b)
{
	GFq_BivariateMonomial result(a);
	result /= b;
	return result;
}

// ================================================================================================
GFq_BivariateMonomial operator /(const GFq_BivariateMonomial& a, const GFq_Element& b)
{
	GFq_BivariateMonomial result(a);
	result /= b;
	return result;
}

// ================================================================================================
GFq_WeightedRevLex_BivariateMonomial::GFq_WeightedRevLex_BivariateMonomial(unsigned int w_x, unsigned int w_y) :
		weights(w_x, w_y)
{}

// ================================================================================================
GFq_WeightedRevLex_BivariateMonomial::GFq_WeightedRevLex_BivariateMonomial(const std::pair<unsigned int, unsigned int>& _weights) :
		weights(_weights)
{}

// ================================================================================================
GFq_WeightedRevLex_BivariateMonomial::~GFq_WeightedRevLex_BivariateMonomial()
{}

// ================================================================================================
bool GFq_WeightedRevLex_BivariateMonomial::operator()(const GFq_BivariateMonomialExponents& e1, const GFq_BivariateMonomialExponents& e2) const
{
	unsigned int w1 = e1.wdeg(weights.first, weights.second);
	unsigned int w2 = e2.wdeg(weights.first, weights.second);

	if (w1 == w2)
	{
		return e1.x() > e2.x();
	}
	else
	{
		return (w1 < w2);
	}
}

// ================================================================================================
std::ostream& operator <<(std::ostream& os, const GFq_BivariateMonomialKeyValueRepresentation& monomial)
{

	if (monomial.first.are_zero())
	{
		os << monomial.second;
	}
	else
	{
		if (monomial.second != 1)
		{
			os << monomial.second;
			os << "*";
		}
	}

	if (monomial.first.x() > 0)
	{
		os << "X";

		if (monomial.first.x() > 1)
		{
			os << "^" << monomial.first.x();
		}

		if (monomial.first.y() > 0)
		{
			os << "*";
		}
	}

	if (monomial.first.y() > 0)
	{
		os << "Y";

		if (monomial.first.y() > 1)
		{
			os << "^" << monomial.first.y();
		}
	}

	return os;
}

// ================================================================================================
GFq_BivariateMonomialKeyValueRepresentation make_bivariate_monomial(GFq_Element coeff, unsigned int exp_x, unsigned int exp_y)
{
	return std::make_pair(std::make_pair(exp_x, exp_y), coeff);
}


} // namespace gf
} // namsepace rssoft
