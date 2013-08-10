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

	 Univariate polynomials with coefficients in GF(2) class hence
	 members of GF(2)[X]

*/

#ifndef __GF2_POLYNOMIAL_H__
#define __GF2_POLYNOMIAL_H__

#include "GF2_Element.h"
#include <iostream>
#include <vector>
#include <utility>
#include <assert.h>

namespace rssoft
{
namespace gf
{

class GF2_Polynomial
{
public:
	/**
	 * Default constructor. Makes an empty thus invalid polynomial
	 */
	GF2_Polynomial();

	/**
	 * Constructs a new polynomial with coefficients in GF(2) with the specified coefficients
	 * \param size The number of coefficients
	 * \param gfe Coefficients array in increasing powers of variable. If null it assign all zero values to coefficients
	 */
	GF2_Polynomial(unsigned int size, GF2_Element* gfe);

	/**
	 * Copy constructor
	 * \param polynomial Reference of the polynomial to copy
	 */
	GF2_Polynomial(const GF2_Polynomial& polynomial);

	/**
	 * Constructs a zero degree polynomial with the given element as coefficient.
	 */
	GF2_Polynomial(const GF2_Element& gfe);

	/**
	 * Constructs the X^n polynomial
	 */
	GF2_Polynomial(unsigned int n);

	/**
	 * Destructor
	 */
	~GF2_Polynomial()
	{};

	/**
	 * Tells if the polynomial has coefficients
	 */
	bool valid() const;

	/**
	 * Tells if the polynomial has no coefficients or a null coefficient
	 */
	bool null() const;

	/**
	 * Degree of the polynomial
	 */
	unsigned int deg() const;

	/**
	 * Gets a constant reference to the coefficients vector
	 */
	const std::vector<GF2_Element>& get_poly() const;

	/**
	 * Gets a modifiable reference to the coefficients vector
	 */
	std::vector<GF2_Element>& get_poly_for_update();

	/**
	 * Sets the degree of the polynomial.
	 * Truncates or extends with zeros the coefficients vector internally
	 */
	void set_degree(unsigned int x);

	GF2_Polynomial& operator =(const GF2_Polynomial& polynomial);
	GF2_Polynomial& operator =(const GF2_Element& gfe);
	GF2_Polynomial& operator+=(const GF2_Polynomial& polynomial);
	GF2_Polynomial& operator+=(const GF2_Element& gfe);
	GF2_Polynomial& operator-=(const GF2_Polynomial& polynomial);
	GF2_Polynomial& operator-=(const GF2_Element& gfe);
	GF2_Polynomial& operator*=(const GF2_Polynomial& polynomial);
	GF2_Polynomial& operator*=(const GF2_Element& gfe);
	GF2_Polynomial& operator/=(const GF2_Polynomial& divisor);
	GF2_Polynomial& operator/=(const GF2_Element& gfe);
	GF2_Polynomial& operator%=(const GF2_Polynomial& divisor);
	GF2_Polynomial& operator%=(const unsigned int& power);
	GF2_Polynomial& operator^=(const int& n);
	GF2_Polynomial& operator<<=(const unsigned int& n);
	GF2_Polynomial& operator>>=(const unsigned int& n);

	GF2_Element& operator[](unsigned int term);
	GF2_Element operator()(GF2_Element value);

	const GF2_Element& operator[](unsigned int term) const;
	const GF2_Element operator()(GF2_Element value) const;

	bool operator==(const GF2_Polynomial& polynomial) const;
	bool operator!=(const GF2_Polynomial& polynomial) const;

	/**
	 * Prints a polynomial to an output stream
	 */
	friend std::ostream& operator <<(std::ostream& os, const GF2_Polynomial& polynomial);

private:
	std::vector<GF2_Element> poly; //!< Coefficients vector
};

void simplify(GF2_Polynomial& polynomial);
GF2_Polynomial operator +(const GF2_Polynomial& a, const GF2_Polynomial& b);
GF2_Polynomial operator +(const GF2_Polynomial& a, const GF2_Element& b);
GF2_Polynomial operator +(const GF2_Element& a, const GF2_Polynomial& b);
GF2_Polynomial operator -(const GF2_Polynomial& a, const GF2_Polynomial& b);
GF2_Polynomial operator -(const GF2_Polynomial& a, const GF2_Element& b);
GF2_Polynomial operator -(const GF2_Element& a, const GF2_Polynomial& b);
GF2_Polynomial operator *(const GF2_Polynomial& a, const GF2_Polynomial& b);
GF2_Polynomial operator *(const GF2_Element& a, const GF2_Polynomial& b);
GF2_Polynomial operator *(const GF2_Polynomial& a, const GF2_Element& b);
GF2_Polynomial operator /(const GF2_Polynomial& a, const GF2_Polynomial& b);
GF2_Polynomial operator /(const GF2_Polynomial& a, const GF2_Element& b);
GF2_Polynomial operator %(const GF2_Polynomial& a, const GF2_Polynomial& b);
GF2_Polynomial operator %(const GF2_Polynomial& a, const unsigned int& power);
GF2_Polynomial operator ^(const GF2_Polynomial& a, const int& n);

/**
 * Shifts coefficients up by n. This increases the degree by n. Lower positions are filled with zero coefficients.
 */
GF2_Polynomial operator <<(const GF2_Polynomial& a, const unsigned int& n);

/**
 * Shifts coefficients up down n. This decreases the degree by n. if n is larger than the polynomial degree this
 * will effectively create an invalid polynomial (without coefficients).
 */
GF2_Polynomial operator >>(const GF2_Polynomial& a, const unsigned int& n);

/**
 * Get the Greatest Common Divisor of two polynomials. The order in which you enter them does not matter. It will
 * make appropriate choices so gcd(a,b) and gcd(b,a) are equivalent. gcd(0,0) is invalid, gcd(a,0) yields a and gcd(0,b)
 * yields b.
 * \param a Polynomial a
 * \param b Polynomial b
 * \return gcd(a,b)
 */
GF2_Polynomial gcd(const GF2_Polynomial& a, const GF2_Polynomial& b);

/**
 * Long (or Euclidean) division returning both quotient and remainder
 * \param a Dividend
 * \param b Divisor
 * \return A pair of polynomials which first member is the quotient and second member is the remainder
 */
std::pair<GF2_Polynomial, GF2_Polynomial> div(const GF2_Polynomial& a,
		const GF2_Polynomial& b);

/**
 * Tells if the polynomial is irreducible
 * \param a Polynomial to test
 * \return true if irreducible
 */
bool irreducible(const GF2_Polynomial& a);

/**
 * Gives the number of non null coefficients of a polynomial
 * \param a Polynomial
 * \return Number of non null coefficients
 */
unsigned int coeff_parity(const GF2_Polynomial& a);

/**
 * Tells if the polynomial is primitive in GF(2^m). Uses heuristics and can serve to eliminate most
 * non primitive polynomials however to build a GF(2^m) field it is recommended to use well known
 * primitive polynomials
 * \param a Polynomial to test
 * \param m m as in GF(2^m)
 * \return true if primitive
 */
bool primitive(const GF2_Polynomial& a, unsigned int m);


} // namespace gf
} // namespace rssoft

#endif // __GF2_POLYNOMIAL_H__

