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
#ifndef __RS_SYSTEMATIC_ENCODING_H__
#define __RS_SYSTEMATIC_ENCODING_H__

#include "GFq.h"
#include "GFq_Element.h"
#include "GFq_Polynomial.h"
#include <vector>

namespace rssoft
{

/**
 * \brief Does the Reed-Solomon systematic encoding of a message. his is the genuine, straightforward, non-systematic encoding that
 * takes message symbols to build the successive coefficients of the encoding polynomial. Then this polynomial is evaluated at the
 * evaluation points to make the codeword.
 */
class RS_SystematicEncoding
{
public:
	/**
	 * Constructor
	 * \param _gf Galois Field in use
	 * \param _k k as in RS(n,k). n is the "size" of the Galois Field
	 * \param _init_power Initial power of alpha
	 */
	RS_SystematicEncoding(const gf::GFq& _gf, unsigned int _k, unsigned int _init_power);

	/**
	 * Destructor. Nothing special.
	 */
	~RS_SystematicEncoding();

	/**
	 * Runs an encoding
	 * \param message Message symbols to be encoded
	 * \param codeword RS codeword in elements of GFq that will be built
	 */
	void run(const std::vector<gf::GFq_Symbol>& message, std::vector<gf::GFq_Symbol>& codeword) const;

protected:
	const gf::GFq& gf; //!< Galois Field in use
	unsigned int k; //!< k as in RS(n,k). n is the "size" of the Galois Field
	unsigned int init_power; //!< Initial power of alpha
	gf::GFq_Polynomial G; //!< Generator polynomial
};


} // namespace rssoft

#endif // __RS_SYSTEMATIC_ENCODING_H__

