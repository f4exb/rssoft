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

#include "EvaluationValues.h"
#include "GFq.h"
#include "RSSoft_Exception.h"

namespace rssoft
{

// ================================================================================================
EvaluationValues::EvaluationValues(const gf::GFq& _gf) :
	gf(_gf)
{
	// default X interpolation values initialized as increasing powers of alpha starting at a^0 = 1
	// default Y interpolation values initialized as the increasing natural order of symbols

	y_values.push_back(gf::GFq_Element(gf, 0));

	for (unsigned int i=0; i < gf.size(); i++)
	{
		x_values.push_back(gf::GFq_Element(gf, gf.alpha(i)));
		y_values.push_back(gf::GFq_Element(gf, i+1));
	}
}

// ================================================================================================
EvaluationValues::EvaluationValues(const gf::GFq& _gf, const std::vector<gf::GFq_Element>& _x_values, const std::vector<gf::GFq_Element>& _y_values) :
	gf(_gf),
	x_values(_x_values),
	y_values(_y_values)
{
	if (x_values.size() > gf.size())
	{
		throw RSSoft_Exception("number of evaluation points cannot be more than the number of non null elements in the field");
	}
	else if (y_values.size() > gf.size()+1)
	{
		throw RSSoft_Exception("number of symbols cannot be more than the number of elements in the field");
	}
}

// ================================================================================================
EvaluationValues::~EvaluationValues()
{}

} // namespace rssoft
