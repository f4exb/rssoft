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

	 Evaluation points of encoding polynomial and symbol values
	 Following the format of Reliability and Multiplicity matrices:
	 - Horizontal, column or X values represent the successive evaluation points in GFq of the
	   encoding polynomial
	 - Vertical, row or Y values represent the successive symbol values of the corresponding elements in GFq

*/
#ifndef __EVALUATION_VALUES__
#define __EVALUATION_VALUES__

#include "GFq_Element.h"
#include <vector>

namespace rssoft
{
namespace gf
{
class GFq;
}


/**
 * \brief Evaluation points of encoding polynomial and symbol values. Following the format of Reliability and Multiplicity matrices:
 * - Horizontal, column or X values represent the successive evaluation points in GFq of the encoding polynomial
 * - Vertical, row or Y values represent the successive symbol values of the corresponding elements in GFq
 */
class EvaluationValues
{
public:
	/**
	 * Default constructor
	 * X values are the successive powers of alpha
	 * Y values are 0 followed by the successive powers of alpha
	 * \param gf Reference to the Galois Field being used
	 */
	EvaluationValues(const gf::GFq& _gf);

	/**
	 * Constructor given X and Y values. It is programmer's responsibility to enter values correctly
	 * \param gf Reference to the Galois Field being used
	 * \param _x_values Successive evaluation points in GFq of the encoding polynomial
	 * \param _y_values Successive symbol values of the corresponding elements in GFq
	 */
	EvaluationValues(const gf::GFq& _gf, const std::vector<gf::GFq_Element>& _x_values, const std::vector<gf::GFq_Element>& _y_values);

	/**
	 * Destructor
	 */
	~EvaluationValues();

	const std::vector<gf::GFq_Element>& get_x_values() const
	{
		return x_values;
	}

	const std::vector<gf::GFq_Element>& get_evaluation_points() const
	{
		return x_values;
	}

	const std::vector<gf::GFq_Element>& get_y_values() const
	{
		return y_values;
	}

	const std::vector<gf::GFq_Element>& get_symbols() const
	{
		return y_values;
	}


protected:
	const gf::GFq& gf; //!< Galois Field being used
	std::vector<gf::GFq_Element> x_values; //!< successive evaluation points in GFq of the encoding polynomial
	std::vector<gf::GFq_Element> y_values; //!< successive symbol values of the corresponding elements in GFq
};

} // namespace rssoft

#endif // __EVALUATION_VALUES__
