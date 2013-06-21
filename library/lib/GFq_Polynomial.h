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
	 Extensively modified augmented with a few methods and included in
	 RSSoft in gf sub-namespace.

	 Univariate polynomials with coefficients in GF(2^m) class hence
	 members of GF(2^m)[X]

*/
/*
 *********************************************************************
 *                                                                   *
 *        Galois Field Arithmetic Library (version 0.0.1)            *
 *                                                                   *
 * Class: Galois Field Polynomial                                    *
 * Version: 0.0.1                                                    *
 * Author: Arash Partow - 2000                                       *
 * URL: http://www.partow.net/projects/galois/index.html             *
 *                                                                   *
 * Copyright Notice:                                                 *
 * Free use of this library is permitted under the guidelines and    *
 * in accordance with the most current version of the Common Public  *
 * License.                                                          *
 * http://www.opensource.org/licenses/cpl.php                        *
 *                                                                   *
 *********************************************************************
 */

#ifndef __GFQ_POLYNOMIAL_H__
#define __GFQ_POLYNOMIAL_H__

#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <assert.h>
#include "GFq.h"
#include "GFq_Element.h"

namespace rssoft
{
namespace gf
{

/**
 * \brief Univariate polynomials with coefficients in GF(2^m)
 */
class GFq_Polynomial
{

public:

	/**
	 * Constructs a new polynomial with coefficients in the specified GF(2^m)
	 * \param _gf The coefficients GF(2^m)
	 */
	GFq_Polynomial(const GFq& _gf);

	/**
	 * Constructs a new polynomial with coefficients in the specified GF(2^m) with the specified coefficients
	 * \param _gf The coefficients GF(2^m)
	 * \param size The number of coefficients
	 * \param gfe Coefficients array in increasing powers of variable
	 */
	GFq_Polynomial(const GFq& _gf, const unsigned int size, GFq_Element* gfe = NULL);

	/**
	 * Constructs a new polynomial with coefficients in the specified GF(2^m) with the specified coefficients
	 * \param _gf The coefficients GF(2^m)
	 * \param gfe Coefficients vector in increasing powers of variable
	 */
	GFq_Polynomial(const GFq& _gf, const std::vector<GFq_Element>& gfe);

	/**
	 * Copy constructor
	 * \param polynomial Reference of the polynomial to copy
	 */
	GFq_Polynomial(const GFq_Polynomial& polynomial);

	/**
	 * Constructs a zero degree polynomial with the given element as coefficient. The reference Galois Field is taken
	 * from the element
	 */
	GFq_Polynomial(const GFq_Element& gfe);

	/**
	 * Destructor
	 */
	~GFq_Polynomial()
	{
	}
	;

	/**
	 * Initializes polynomial with a vector of coefficients in increasing power of variable
	 */
	void init(const std::vector<GFq_Element>& _poly);

	/**
	 * Tells if the polynomial has coefficients
	 */
	bool is_valid() const;

	/**
	 * Tells if the polynomial has no coefficients or a null coefficient P(X)=0
	 */
	bool is_zero() const;

	/**
	 * Tells if the polynomial is P(X)=1
	 */
	bool is_one() const;

	/**
	 * Degree of the polynomial
	 */
	unsigned int deg() const;

	/**
	 * Reference to the Galois Field of the coefficients of the polynomial
	 */
	const GFq& field() const;

	/**
	 * Gets a constant reference to the coefficients vector
	 */
	const std::vector<GFq_Element>& get_poly() const;

	/**
	 * Gets a modifiable reference to the coefficients vector
	 */
	std::vector<GFq_Element>& get_poly_for_update();

	/**
	 * Gets the coefficients vector as symbols
	 */
    void get_poly_symbols(std::vector<GFq_Symbol>& symbols, unsigned int size=0) const
    {
        symbols.resize((size == 0 ? poly.size() : size));
        std::transform(poly.begin(), poly.end(), symbols.begin(), gfq_element_to_symbol);
    }

	/**
	 * Sets the degree of the polynomial.
	 * Truncates or extends with zeros the coefficients vector internally
	 */
	void set_degree(const unsigned int& x);

	GFq_Polynomial& operator =(const GFq_Polynomial& polynomial);
	GFq_Polynomial& operator =(const GFq_Element& gfe);
	GFq_Polynomial& operator+=(const GFq_Polynomial& polynomial);
	GFq_Polynomial& operator+=(const GFq_Element& gfe);
	GFq_Polynomial& operator-=(const GFq_Polynomial& polynomial);
	GFq_Polynomial& operator-=(const GFq_Element& gfe);
	GFq_Polynomial& operator*=(const GFq_Polynomial& polynomial);
	GFq_Polynomial& operator*=(const GFq_Element& gfe);
	GFq_Polynomial& operator/=(const GFq_Polynomial& divisor);
	GFq_Polynomial& operator/=(const GFq_Element& gfe);
	GFq_Polynomial& operator%=(const GFq_Polynomial& divisor);
	GFq_Polynomial& operator%=(const unsigned int& power);
	GFq_Polynomial& operator^=(const int& n);
	GFq_Polynomial& operator<<=(const unsigned int& n);
	GFq_Polynomial& operator>>=(const unsigned int& n);

	GFq_Element& operator[](const unsigned int& term);
	GFq_Element operator()(const GFq_Element& value);
	GFq_Element operator()(GFq_Symbol value);

	const GFq_Element& operator[](const unsigned int& term) const;
	const GFq_Element operator()(const GFq_Element& value) const;
	const GFq_Element operator()(GFq_Symbol value) const;

	bool operator==(const GFq_Polynomial& polynomial) const;
	bool operator!=(const GFq_Polynomial& polynomial) const;

	/**
	 * Calculates the derivative of a polynomial
	 */
	GFq_Polynomial derivative();
    
    /**
     * Make this polynomial monic
     * \return the leading coefficient by which all coefficients are divided
     */
    GFq_Element make_monic();

    /**
     * Find non null roots of polynomial by Chien search. Basically this is an exhaustive search
     * optimized for hardware implementation but is also interesting in software.
     * Includes null root afterwards if constant coefficient is zero.
     * see: http://www.stanford.edu/class/ee387/handouts/notes7.pdf
     * \param roots Vector of root field elements filled by the method
     */
     void rootChien(std::vector<GFq_Element>& roots);

	/**
	 * Sets the output of coefficients in a power of alpha (printed a^i) ot binary representation
	 * \param _alpha_format true: power of alpha representation. false: binary representation
	 */
	void set_alpha_format(bool _alpha_format)
	{
		alpha_format = _alpha_format;
	}

	/**
	 * Prints a polynomial to an output stream
	 */
	friend std::ostream& operator <<(std::ostream& os, const GFq_Polynomial& polynomial);

protected:
	const GFq& gf; //!< Reference to the Galois Field of the coefficients
	std::vector<GFq_Element> poly; //!< Coefficients vector
	bool alpha_format; // true: power of alpha representation, false: binary representation of printed coefficients
};

void simplify(GFq_Polynomial& polynomial);
GFq_Polynomial operator +(const GFq_Polynomial& a, const GFq_Polynomial& b);
GFq_Polynomial operator +(const GFq_Polynomial& a, const GFq_Element& b);
GFq_Polynomial operator +(const GFq_Element& a, const GFq_Polynomial& b);
GFq_Polynomial operator +(const GFq_Polynomial& a, const GFq_Symbol& b);
GFq_Polynomial operator +(const GFq_Symbol& a, const GFq_Polynomial& b);
GFq_Polynomial operator -(const GFq_Polynomial& a, const GFq_Polynomial& b);
GFq_Polynomial operator -(const GFq_Polynomial& a, const GFq_Element& b);
GFq_Polynomial operator -(const GFq_Element& a, const GFq_Polynomial& b);
GFq_Polynomial operator -(const GFq_Polynomial& a, const GFq_Symbol& b);
GFq_Polynomial operator -(const GFq_Symbol& a, const GFq_Polynomial& b);
GFq_Polynomial operator *(const GFq_Polynomial& a, const GFq_Polynomial& b);
GFq_Polynomial operator *(const GFq_Element& a, const GFq_Polynomial& b);
GFq_Polynomial operator *(const GFq_Polynomial& a, const GFq_Element& b);
GFq_Polynomial operator /(const GFq_Polynomial& a, const GFq_Polynomial& b);
GFq_Polynomial operator /(const GFq_Polynomial& a, const GFq_Element& b);
GFq_Polynomial operator %(const GFq_Polynomial& a, const GFq_Polynomial& b);
GFq_Polynomial operator %(const GFq_Polynomial& a, const unsigned int& power);
GFq_Polynomial operator ^(const GFq_Polynomial& a, const int& n);

/**
 * Shifts coefficients up by n. This increases the degree by n. Lower positions are filled with zero coefficients.
 */
GFq_Polynomial operator <<(const GFq_Polynomial& a, const unsigned int& n);

/**
 * Shifts coefficients up down n. This decreases the degree by n. if n is larger than the polynomial degree this
 * will effectively create an invalid polynomial (without coefficients).
 */
GFq_Polynomial operator >>(const GFq_Polynomial& a, const unsigned int& n);

/**
 * Get the Greatest Common Divisor of two polynomials. The order in which you enter them does not matter. It will
 * make appropriate choices so gcd(a,b) and gcd(b,a) are equivalent. gcd(0,0) is invalid, gcd(a,0) yields a and gcd(0,b)
 * yields b.
 * \param a Polynomial a
 * \param b Polynomial b
 * \return gcd(a,b)
 */
GFq_Polynomial gcd(const GFq_Polynomial& a, const GFq_Polynomial& b);

/**
 * Long (or Euclidean) division returning both quotient and remainder
 * \return A pair of polynomials which first member is the quotient and second member is the remainder
 */
std::pair<GFq_Polynomial, GFq_Polynomial> div(const GFq_Polynomial& a,
		const GFq_Polynomial& b);

/**
 * Find non null roots of polynomial by exhaustive search
 * \param a Polynomial which roots are searched
 * \return vector of root field elements
 */
std::vector<GFq_Element> rootex_nz(const GFq_Polynomial& a);

/**
 * Find roots of polynomial by exhaustive search
 * \param a Polynomial which roots are searched
 * \return vector of root field elements
 */
std::vector<GFq_Element> rootex(const GFq_Polynomial& a);

/** 
 * Get the monic polynomial from a polynomial
 * \param a Polynomial to get monic from
 * \param lead_coeff Reference to a GaloisFieldElement which will receive the original leading polynomial
 * \return The monic polynomial
 */
GFq_Polynomial get_monic(const GFq_Polynomial& a, GFq_Element& lead_poly);

/**
 * Square free decomposition of a polynomial using Yun's algorithm
 * WARNING: largely untested and not validated
 * \param a Polynomial from which to do the square free decomposition
 * \return the vector of pairs of polynomial factors and their exponent
 */
std::vector<GFq_Polynomial> square_free_decomposition(const GFq_Polynomial& a);

} // namespace gf
} // namespace rssoft

#endif // __GFQ_POLYNOMIAL_H__
