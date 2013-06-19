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

 Reed-Solomon encoding class

 */
#include "RS_Encoding.h"
#include "GFq.h"
#include "GFq_Polynomial.h"
#include "EvaluationValues.h"
#include "RSSoft_Exception.h"

namespace rssoft
{

// ================================================================================================
RS_Encoding::RS_Encoding(const gf::GFq& _gf, unsigned int _k, const EvaluationValues& _evaluation_values) :
	gf(_gf),
	k(_k),
	evaluation_values(_evaluation_values)
{}

// ================================================================================================
RS_Encoding::~RS_Encoding()
{}

// ================================================================================================
void RS_Encoding::run(const std::vector<gf::GFq_Symbol>& message, std::vector<gf::GFq_Symbol>& codeword) const
{
	if (message.size() != k)
	{
		throw RSSoft_Exception("Invalid message length");
	}
	else
	{
		std::vector<gf::GFq_Symbol>::const_iterator s_it = message.begin();
		std::vector<gf::GFq_Element> encoding_coefficients;

		for(; s_it != message.end(); ++s_it)
		{
			encoding_coefficients.push_back(gf::GFq_Element(gf, *s_it));
		}

		gf::GFq_Polynomial encoding_polynomial(gf, encoding_coefficients);
		codeword.clear();
		const std::vector<gf::GFq_Element>& evaluation_points = evaluation_values.get_evaluation_points();
		std::vector<gf::GFq_Element>::const_iterator evp_it = evaluation_points.begin();

		for (; evp_it != evaluation_points.end(); ++evp_it)
		{
			codeword.push_back(encoding_polynomial(*evp_it).poly());
		}
	}
}

} // namespace rssoft
