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

 Guruswami-Sudan-Koetter-Vardy interpolation class for soft decision
 decoding

 */
#ifndef __GSKV_INTERPOLATION_H__
#define __GSKV_INTERPOLATION_H__

#include "GFq_BivariatePolynomial.h"
#include "GFq_Element.h"
#include <utility>
#include <vector>

namespace rssoft
{

namespace gf
{
class GFq;
}

class MultiplicityMatrix;

class GSKV_Interpolation
{
public:
	/**
	 * Constructor
	 * \param gf Reference to the Galois Field being used
	 * \param k as in RS(n,k)
	 */
	GSKV_Interpolation(const gf::GFq& _gf, unsigned int _k);

	/**
	 * Destructor
	 */
	~GSKV_Interpolation();

	/**
	 * Run the interpolation based on given multiplicity matrix
	 */
	void run(const MultiplicityMatrix& mmat);

protected:
	/**
	 * Interpolation polynomial maximum degrees
	 * \return (dX,dY) pair of maximum degrees in X and Y
	 */
	std::pair<unsigned int, unsigned int> maximum_degrees(const MultiplicityMatrix& mmat);

	/**
	 * Initialize G list of polynomials and related lists
	 */
	void init_G(unsigned int dY);

	/**
	 * Process an interpolation point with multiplicity. This is the outer iteration of the algorithm
	 * \param iX X coordinate that is evaluation point in GFq
	 * \param iY Y coordinate that is value in GFq at evaluation point
	 * \param multiplicity Multiplicity
	 */
	void process_point(unsigned int iX, unsigned int iY, unsigned int multiplicity);

	/**
	 * Process a Hasse derivative. This is the inner iteration of the algorithm
	 * \param x Evaluation point in GFq
	 * \param y Value in GFq at evaluation point
	 * \param mu Mu parameter of Hasse derivative (related to X)
	 * \param nu Nu parameter of Hasse derivative (related to Y)
	 */
	void process_hasse(const gf::GFq_Element& x, const gf::GFq_Element& y, unsigned int mu, unsigned int nu);

	/**
	 * Finalize process with G list of polynomials and find result polynomial
	 * \return Reference to the result polynomial
	 */
	const gf::GFq_BivariatePolynomial& final_G();

	// fixed parameters
	const gf::GFq& gf; //!< Reference to the Galois Field being used
	unsigned int k; //!< k factor as in RS(n,k)
	std::vector<gf::GFq_Element> x_values; //!< Interpolation X values
	std::vector<gf::GFq_Element> y_values; //!< Interpolation Y values

	// parameters changing at each process run
	std::vector<gf::GFq_BivariatePolynomial> G; //!< The G list of polynomials
	std::vector<bool> calcG; //!< Li Chen's optimization. If true the corresponding polynomial in G is processed.
	std::vector<unsigned int> lodG; //!< Leading orders of polynomials in G
    unsigned int it_number; //!< Hasse derivative iteration number (inner loop)
    unsigned int Cm; //!< Cost of current multiplicity matrix
    unsigned int final_ig; //!< Index of the result polynomial in G list
};

} // namespace rssoft

#endif // __GSKV_INTERPOLATION_H__
