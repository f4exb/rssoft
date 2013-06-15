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
#include "GFq_Polynomial.h"
#include <map>

namespace rssoft
{
namespace gf
{

/**
 * \brief Bivariate polynomials with coefficients in GF(2^m) thus the polynomials in
 * GF(2^m)[X,Y]. Division is intentionally omitted as this is complex and unnecessary for
 * the application
 */
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
	 * Constructs a new empty (thus invalid) bivariate polynomial
	 * \param weights Weight in X,Y pair for monomials weighted ordering
	 */
	GFq_BivariatePolynomial(const std::pair<unsigned int, unsigned int>& _weights);

	/**
	 * Copy constructor.
	 * \param polynomial Polynomial to copy from
	 */
	GFq_BivariatePolynomial(const GFq_BivariatePolynomial& polynomial);

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
	 * Initializes the polynomial with the monomials of another polynomial
	 */
	void init(const GFq_BivariatePolynomial& polynomial);

	/**
	 * Initializes the polynomial with a map of monomials
	 */
	void init(const std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& monomials);
    
    /**
     * Initializes the polynomial as X^n
     * \param gf Galois Field to create the unit GF element as coefficient
     * \param max_pow power of polynomial's unique monomial
     */
    void init_x_pow(const GFq& gf, unsigned int x_pow);

    /**
     * Initializes the polynomial as Y^n
     * \param gf Galois Field to create the unit GF element as coefficient
     * \param max_pow power of polynomial's unique monomial
     */
    void init_y_pow(const GFq& gf, unsigned int y_pow);

    /**
     * Initializes the polynomial as 1+X+..+X^n
     * \param gf Galois Field to create the unit GF element as coefficient
     * \param max_pow Maximum power in the series
     */
    void init_x_pow_series(const GFq& gf, unsigned int max_pow);

    /**
     * Initializes the polynomial as 1+Y+..+Y^n
     * \param gf Galois Field to create the unit GF element as coefficient
     * \param max_pow Maximum power in the series
     */
    void init_y_pow_series(const GFq& gf, unsigned int max_pow);

	/**
	 * Gets the weights pair in (X,Y) used for weighted monomial ordering
	 * \return Reference to the weights pair in (X,Y) used for weighted monomial ordering
	 */
	const std::pair<unsigned int, unsigned int>& get_weights() const
	{
		return weights;
	}

	/**
	 * Tells if the polynomial has coefficients
	 */
	bool is_valid() const;

	/**
	 * Tells if the polynomial has only one constant coefficient with the given value
	 */
	bool is_const(GFq_Element& const_value) const;

	/**
	 * Tells if the polynomial has no coefficients or a null coefficient
	 */
	bool is_zero() const;

	/**
	 * Tells if the polynomial is P(X)=1
	 */
	bool is_one() const;

	/**
	 * Get the monomials
	 * \return Reference to the set of monomials
	 */
	const std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& get_monomials() const
	{
		return monomials;
	}

	/**
	 * Get the monomials to be updated
	 * \return Reference to the set of monomials
	 */
	std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& get_monomials_for_update() 
	{
		return monomials;
	}
    
	/**
	 * Get a copy of the leading monomial with respect to weighted reverse lexical order
	 * \return Reference to the leading monomial
	 */
	GFq_BivariateMonomial get_leading_monomial() const
	{
		return static_cast<GFq_BivariateMonomial>(*(monomials.rbegin()));
	}
    
    /**
     * Get X power of leading monomial
     */
    unsigned int lmX() const
    {
        return static_cast<GFq_BivariateMonomial>(*(monomials.rbegin())).eX();
    }
    
    /**
     * Get Y power of leading monomial
     */
    unsigned int lmY() const
    {
        return static_cast<GFq_BivariateMonomial>(*(monomials.rbegin())).eY();
    }
    
    /**
     * Weighted degree of polynomial. That is the highest weighted degree of its monomials. 
     */
    unsigned int wdeg() const;

	/**
	 * Helper method to create the vector of monomials of the sum of polynomials a and b
	 */
	static void sum(std::vector<GFq_BivariateMonomial>& sum_monomials, const GFq_BivariatePolynomial& a, const GFq_BivariatePolynomial& b);

	/**
	 * Helper method to create the map of monomials of the product of polynomials a and b
	 */
	static void product(std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& prod_monomials,
			const GFq_BivariatePolynomial& a,
			const GFq_BivariatePolynomial& b);

	/**
	 * Helper method to create the map of monomials of the division of polynomial a by monomial b
	 */
	static void division(std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& div_monomials,
			const GFq_BivariatePolynomial& a,
			const GFq_BivariateMonomial& b);

	GFq_BivariatePolynomial& operator =(const GFq_BivariatePolynomial& polynomial);
	GFq_BivariatePolynomial& operator+=(const GFq_BivariatePolynomial& polynomial);
	GFq_BivariatePolynomial& operator+=(const GFq_Element& gfe);
	GFq_BivariatePolynomial& operator-=(const GFq_BivariatePolynomial& polynomial);
	GFq_BivariatePolynomial& operator-=(const GFq_Element& gfe);
	GFq_BivariatePolynomial& operator*=(const GFq_BivariatePolynomial& polynomial);
	GFq_BivariatePolynomial& operator*=(const GFq_BivariateMonomial& monomial);
	GFq_BivariatePolynomial& operator*=(const GFq_Element& gfe);
	GFq_BivariatePolynomial& operator/=(const GFq_BivariateMonomial& monomial);
	GFq_BivariatePolynomial& operator/=(const GFq_Element& gfe);

	bool operator==(const GFq_BivariatePolynomial& polynomial) const;
	bool operator!=(const GFq_BivariatePolynomial& polynomial) const;

	/**
	 * Evaluation of bivariate polynomial at a (x,y) point in GFq^2
	 * \param x_value x coordinate
	 * \param y_value y coordinate
	 * \return Value of polynomial at (x,y) point
	 */
	GFq_Element operator()(const GFq_Element& x_value, const GFq_Element& y_value) const;

	/**
	 * Evaluation of polynomial for Y=0 as a univariate polynomial in X
	 * \return Univariate polynomial in X as P(X,0)
	 */
	GFq_Polynomial get_X_0() const;
	
	/**
	 * Evaluation of polynomial for X=0 as a univariate polynomial in Y
	 * \return Univariate polynomial in Y as P(0,Y)
	 */
	GFq_Polynomial get_0_Y() const;

	/**
	 * Applies to self the star function as P*(X,Y) = P(X,Y)/X^h where h is the greatest power of X so that X^h divides P
	 * \return reference to the new polynomial
	 */
	GFq_BivariatePolynomial& make_star();

 	/**
 	 * Applies to self the [mu,nu] Hasse derivative
 	 * \param mu mu parameter (applies to X factors)
 	 * \param nu nu parameter (applies to Y factors)
 	 * \return reference to the new polynomial
 	 */
 	GFq_BivariatePolynomial& make_dHasse(unsigned int mu, unsigned int nu);

	/**
	 * Prints a polynomial to an output stream
	 */
	friend std::ostream& operator <<(std::ostream& os, const GFq_BivariatePolynomial& polynomial);

protected:

	GFq_Polynomial get_v_0(bool x_terms) const;
    
	std::pair<unsigned int, unsigned int> weights; //<! weights for weighted degree ordering
	std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial> monomials; //<! set of monomials
};

/**
 * Removes monomials with zero coefficients
 */
void simplify(GFq_BivariatePolynomial& polynomial);

GFq_BivariatePolynomial operator +(const GFq_BivariatePolynomial& a, const GFq_BivariatePolynomial& b);
GFq_BivariatePolynomial operator +(const GFq_BivariatePolynomial& a, const GFq_Element& b);
GFq_BivariatePolynomial operator +(const GFq_Element& a, const GFq_BivariatePolynomial& b);
GFq_BivariatePolynomial operator -(const GFq_BivariatePolynomial& a, const GFq_BivariatePolynomial& b);
GFq_BivariatePolynomial operator -(const GFq_BivariatePolynomial& a, const GFq_Element& b);
GFq_BivariatePolynomial operator -(const GFq_Element& a, const GFq_BivariatePolynomial& b);
GFq_BivariatePolynomial operator *(const GFq_BivariatePolynomial& a, const GFq_BivariatePolynomial& b);
GFq_BivariatePolynomial operator *(const GFq_BivariatePolynomial& a, const GFq_BivariateMonomial& b);
GFq_BivariatePolynomial operator *(const GFq_BivariatePolynomial& a, const GFq_Element& b);
GFq_BivariatePolynomial operator *(const GFq_Element& a, const GFq_BivariatePolynomial& b);
GFq_BivariatePolynomial operator /(const GFq_BivariatePolynomial& a, const GFq_BivariateMonomial& b);
GFq_BivariatePolynomial operator /(const GFq_BivariatePolynomial& a, const GFq_Element& b);

/**
 * Star function as P*(X,Y) = P(X,Y)/X^h where h is the greatest power of X so that X^h divides P
 * \param a Input polynomial
 * \return P*(X,Y)
 */
GFq_BivariatePolynomial star(const GFq_BivariatePolynomial& a);

/**
 * [mu,nu] Hasse derivative
 * \param mu mu parameter (applies to X factors)
 * \param nu nu parameter (applies to Y factors)
 * \param a Input polynomial
 * \return [mu,nu] Hasse derivative of the polynomial
 */
GFq_BivariatePolynomial dHasse(unsigned int mu, unsigned int nu, const GFq_BivariatePolynomial& a);

} // namespace gf
} // namespace rssoft

#endif // __GFQ_BIVARIATE_POLYNOMIAL_H__
