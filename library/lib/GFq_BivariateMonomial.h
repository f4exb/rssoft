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
#ifndef __GFQ_BIVARIATE_MONOMIAL_H__
#define __GFQ_BIVARIATE_MONOMIAL_H__

#include <GFq_Element.h>
#include <utility>

namespace rssoft
{
namespace gf
{

/**
 * \brief Bivariate monomials with coefficient in GF(2^m)
 */
class GFq_BivariateMonomial
{
public:
	/**
	 * Constructs a new bivariate monomial
	 * \param coeff Coefficient
	 * \param x_pow Power in X
	 * \param y_pow Power in Y
	 */
	GFq_BivariateMonomial(const GFq_Element& coeff, unsigned int x_pow, unsigned int y_pow);

	/**
	 * Destructor
	 */
	~GFq_BivariateMonomial();

	/**
	 * Returns power in X
	 */
	unsigned int x_pow() const
	{
		return exponents.first;
	}

	/**
	 * Returns power in Y
	 */
	unsigned int y_pow() const
	{
		return exponents.second;
	}

	/**
	 * Weighted degree of monomial
	 */
	unsigned int wdeg(unsigned int w_x, unsigned int w_y) const
	{
		return w_x*exponents.first + w_y*exponents.second;
	}

	/**
	 * Prints a monomial to an output stream
	 */
	friend std::ostream& operator <<(std::ostream& os, const GFq_BivariateMonomial& monomial);

protected:
	std::pair<unsigned int, unsigned int> exponents;
	GFq_Element coefficient;
};

/**
 * \brief Weighted reverse lexical order of bivariate monomials
 */
class GFq_WeightedRevLex_BivariateMonomial
{
public:
	/**
	 * Constructs a new weighted reverse lexical order of bivariate monomials
	 * \param w_x Weight in X
	 * \param w_y Weight in Y
	 */
	GFq_WeightedRevLex_BivariateMonomial(unsigned int w_x, unsigned int w_y);

	/**
	 * Ordering method
	 * \param m1 First monomial, if lesser than second then return is true
	 * \param m2 Second monomial, if greated than first then return is true
	 */
	bool operator()(const GFq_BivariateMonomial& m1, const GFq_BivariateMonomial& m2) const;

	/**
	 * Destructor
	 */
	~GFq_WeightedRevLex_BivariateMonomial();

protected:
	std::pair<unsigned int, unsigned int> weights;
};

} // namespace gf
} // namespace rssoft

#endif // __GFQ_BIVARIATE_MONOMIAL_H__
