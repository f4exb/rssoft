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
class EvaluationValues;

class GSKV_Interpolation
{
public:
	/**
	 * Constructor
	 * \param _gf Reference to the Galois Field being used
	 * \param _k as in RS(n,k)
	 * \param _evaluation_values Evaluation X,Y values used for coding
	 */
	GSKV_Interpolation(const gf::GFq& _gf, unsigned int _k, const EvaluationValues& _evaluation_values);

	/**
	 * Destructor
	 */
	~GSKV_Interpolation();

	/**
	 * (Re)initialize internal objects to prepare a new run
	 */
	void init();
    
    /**
     * Get the X evaluation points used to calculate candidate codewords (horizontal or column matrix wise)
     */
    const EvaluationValues& get_evaluation_values() const
    {
        return evaluation_values;
    }

    /**
     * Set or reset verbose mode.
     * \param _verbose Verbose level. 0 to shut down any debug message. Active only in debug mode (_DEBUG defined)
     */
    void set_verbosity(unsigned int _verbosity)
    {
        verbosity = _verbosity;
    }
    
    unsigned int get_dX() const
    {
    	return dX;
    }

    unsigned int get_dY() const
    {
    	return dY;
    }

	/**
	 * Run the interpolation based on given multiplicity matrix
     * \return reference to the result polynomial
	 */
	const gf::GFq_BivariatePolynomial& run(const MultiplicityMatrix& mmat);

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
	const EvaluationValues& evaluation_values; //!< Interpolation X,Y values
    unsigned int verbosity; //!< Verbose level, 0 to shut down any debug message

	// parameters changing at each process run
    unsigned int dX;
    unsigned int dY;
    unsigned int mcost; //!< Multiplicity matrix cost
	std::vector<gf::GFq_BivariatePolynomial> G; //!< The G list of polynomials
	std::vector<bool> calcG; //!< Li Chen's optimization. If true the corresponding polynomial in G is processed.
	std::vector<unsigned int> lodG; //!< Leading orders of polynomials in G
    unsigned int it_number; //!< Hasse derivative iteration number (inner loop)
    unsigned int Cm; //!< Cost of current multiplicity matrix
    unsigned int final_ig; //!< Index of the result polynomial in G list
};

} // namespace rssoft

#endif // __GSKV_INTERPOLATION_H__
