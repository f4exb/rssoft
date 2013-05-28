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
#ifndef __GFQ_BIVARIATE_POLYNOMIAL_H__
#define __GFQ_BIVARIATE_POLYNOMIAL_H__

#include "GFq_BivariateMonomial.h"
#include <set>

namespace rssoft
{
namespace gf
{

class GFq_BivariatePolynomial
{
public:
	/**
	 * Constructs a new empty (thus invalid) bivariate polynomial
	 * \param w_x Weight in X for monomials weighted ordering
	 * \param w_y Weight in Y for monomials weighted ordering
	 */
	GFq_BivariatePolynomial(unsigned int w_x, unsigned int w_y);

	/**
	 * Destructor
	 */
	~GFq_BivariatePolynomial();

	/**
	 * Initializes the polynomial as a sum of monomials
	 * \param monomials List of monomials
	 */
	void init(std::vector<GFq_BivariateMonomial>& _monomials);

	/**
	 * Gets the weights pair in (X,Y) used for weighted monomial ordering
	 * \return Reference to the weights pair in (X,Y) used for weighted monomial ordering
	 */
	const std::pair<unsigned int, unsigned int>& get_weights() const
	{
		return weights;
	}

	/**
	 * Gets the leading monomial with respect to ordering
	 * \return Reference to the leading monomial
	 */
	const GFq_BivariateMonomial& leading_monomial() const
	{
		return *(monomials.rbegin());
	}

	/**
	 * Prints a polynomial to an output stream
	 */
	friend std::ostream& operator <<(std::ostream& os, const GFq_BivariatePolynomial& polynomial);


protected:
	std::pair<unsigned int, unsigned int> weights; //<! weights for weighted degree ordering
	std::set<GFq_BivariateMonomial, GFq_WeightedRevLex_BivariateMonomial> monomials; //<! set of monomials
};

} // namespace gf
} // namespace rssoft

#endif // __GFQ_BIVARIATE_POLYNOMIAL_H__
