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
 * \brief Bivariate monomial exponents
 */
class GFq_BivariateMonomialExponents : public std::pair<unsigned int, unsigned int>
{
public:
	GFq_BivariateMonomialExponents(unsigned int eX, unsigned int eY);
	GFq_BivariateMonomialExponents(const std::pair<unsigned int, unsigned int>& _exponents);

	/**
	 * Returns the exponent in x
	 */
	unsigned int x() const
	{
		return first;
	}

	/**
	 * Returns the exponent in y
	 */
	unsigned int y() const
	{
		return second;
	}

	/**
	 * Returns the weighted degree
	 * \param wX Weight in X
	 * \param wY Weight in Y
	 */
	unsigned int wdeg(unsigned int wX, unsigned int wY) const
	{
		return wX*first + wY*second;
	}

	/**
	 * Returns the weighted degree
	 * \param weihgts Pair of weights in (X,Y)
	 */
	unsigned int wdeg(const std::pair<unsigned int, unsigned int>& weights) const
	{
		return weights.first*first + weights.second*second;
	}

	/**
	 * Tests if both exponents are 0 (monomial represents the constant coefficient)
	 */
	bool are_zero() const
	{
		return (first == 0) && (second == 0);
	}
};

/**
 * Monomial representation as a (key, value) pair where key is the pair of exponents and value the coefficient
 * It serves the purpose to be directly usable in the bivariate polynomial's map of monomials
 */
typedef std::pair<GFq_BivariateMonomialExponents, GFq_Element> GFq_BivariateMonomialKeyValueRepresentation;

/**
 * \brief Bivariate Monomial class
 */
class GFq_BivariateMonomial : public std::pair<GFq_BivariateMonomialExponents, GFq_Element>
{
public:
	GFq_BivariateMonomial(const GFq_Element& coeff, unsigned int eX, unsigned int eY);
	GFq_BivariateMonomial(const GFq_Element& coeff, const GFq_BivariateMonomialExponents& exponents);
	GFq_BivariateMonomial(const GFq_BivariateMonomialKeyValueRepresentation& monomial_rep);

	const GFq_Element& coeff() const
	{
		return second;
	}

	unsigned int eX() const
	{
		return first.first;
	}

	unsigned int eY() const
	{
		return first.second;
	}

	const GFq_BivariateMonomialExponents& get_exponents() const
	{
		return first;
	}

	unsigned int wdeg(unsigned int wX, unsigned int wY) const
	{
		return first.wdeg(wX, wY);
	}

	unsigned int wdeg(const std::pair<unsigned int, unsigned int>& weights) const
	{
		return first.wdeg(weights);
	}

	GFq_BivariateMonomial& operator =(const GFq_BivariateMonomial& monomial);
	GFq_BivariateMonomial& operator+=(const GFq_BivariateMonomial& monomial);
	GFq_BivariateMonomial& operator+=(const GFq_Element& gfe);
	GFq_BivariateMonomial& operator-=(const GFq_BivariateMonomial& monomial);
	GFq_BivariateMonomial& operator-=(const GFq_Element& gfe);
	GFq_BivariateMonomial& operator*=(const GFq_BivariateMonomial& monomial);
	GFq_BivariateMonomial& operator*=(const GFq_Element& gfe);
	GFq_BivariateMonomial& operator/=(const GFq_BivariateMonomial& monomial);
	GFq_BivariateMonomial& operator/=(const GFq_Element& gfe);
};

GFq_BivariateMonomial operator +(const GFq_BivariateMonomial& a, const GFq_BivariateMonomial& b);
GFq_BivariateMonomial operator +(const GFq_BivariateMonomial& a, const GFq_Element& b);
GFq_BivariateMonomial operator +(const GFq_Element& a, const GFq_BivariateMonomial& b);
GFq_BivariateMonomial operator -(const GFq_BivariateMonomial& a, const GFq_BivariateMonomial& b);
GFq_BivariateMonomial operator -(const GFq_BivariateMonomial& a, const GFq_Element& b);
GFq_BivariateMonomial operator -(const GFq_Element& a, const GFq_BivariateMonomial& b);
GFq_BivariateMonomial operator *(const GFq_BivariateMonomial& a, const GFq_BivariateMonomial& b);
GFq_BivariateMonomial operator *(const GFq_Element& a, const GFq_BivariateMonomial& b);
GFq_BivariateMonomial operator *(const GFq_BivariateMonomial& a, const GFq_Element& b);
GFq_BivariateMonomial operator /(const GFq_BivariateMonomial& a, const GFq_BivariateMonomial& b);
GFq_BivariateMonomial operator /(const GFq_BivariateMonomial& a, const GFq_Element& b);

/**
 * Prints monomials to output stream
 */
std::ostream& operator <<(std::ostream& os, const GFq_BivariateMonomialKeyValueRepresentation& monomial);

/**
 * Helper to create a monomial representation
 */
GFq_BivariateMonomialKeyValueRepresentation make_bivariate_monomial(GFq_Element coeff, unsigned int exp_x, unsigned int exp_y);

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
	 * Constructs a new weighted reverse lexical order of bivariate monomials
	 * \param weights Weights in (X,Y)
	 */
	GFq_WeightedRevLex_BivariateMonomial(const std::pair<unsigned int, unsigned int>& weights);

	/**
	 * Ordering method
	 * \param e1 First monomial exponents, if lesser than second then return is true
	 * \param e2 Second monomial exponents, if greater than first then return is true
	 */
	bool operator()(const GFq_BivariateMonomialExponents& e1, const GFq_BivariateMonomialExponents& e2) const;

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
