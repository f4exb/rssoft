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

 Does the Reed-Solomon encoding of a message.
 This is the genuine, straightforward, non-systematic encoding that
 takes message symbols to build the successive coefficients of the
 encoding polynomial. Then this polynomial is evaluated at the
 evaluation points to make the codeword.

 */
#ifndef __RS_ENCODING_H__
#define __RS_ENCODING_H__

#include "GFq.h"
#include "GFq_Element.h"
#include <vector>

namespace rssoft
{

class EvaluationValues;

/**
 * \brief Does the Reed-Solomon encoding of a message. his is the genuine, straightforward, non-systematic encoding that
 * takes message symbols to build the successive coefficients of the encoding polynomial. Then this polynomial is evaluated at the
 * evaluation points to make the codeword.
 */
class RS_Encoding
{
public:
	/**
	 * Constructor
	 * \param _gf Galois Field in use
	 * \param _k k as in RS(n,k). n is the "size" of the Galois Field
	 * \param _evaluation_values Evaluation X,Y values of the code
	 */
	RS_Encoding(const gf::GFq& _gf, unsigned int k, const EvaluationValues& _evaluation_values);

	/**
	 * Destructor. Nothing special.
	 */
	~RS_Encoding();

	/**
	 * Runs an encoding
	 * \param message Message symbols to be encoded
	 * \param codeword RS codeword in elements of GFq that will be built
	 */
	void run(const std::vector<gf::GFq_Symbol>& message, std::vector<gf::GFq_Symbol>& codeword) const;

	/**
	 * Runs a systematic encoding
	 * \param initial power of alpha
	 * \param message Message symbols to be encoded
	 * \param codeword RS codeword in elements of GFq that will be built
	 */
	void run_systematic(unsigned int init_pow, const std::vector<gf::GFq_Symbol>& message, std::vector<gf::GFq_Symbol>& codeword) const;

protected:
	const gf::GFq& gf; //!< Galois Field in use
	unsigned int k; //!< k as in RS(n,k). n is the "size" of the Galois Field
	const EvaluationValues& evaluation_values; //!< Evaluation X,Y values of the code
};


} // namespace rssoft

#endif // __RS_ENCODING_H__

