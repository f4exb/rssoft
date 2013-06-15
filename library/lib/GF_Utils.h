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

 Utility routines not linked to a patricular GF entity

 */

#ifndef __GF_UTILS_H__
#define __GF_UTILS_H__

namespace rssoft
{
namespace gf
{

/**
 * Computes parity of a binomial coefficient
 * \return true if binomial coefficient is even
 */
bool binomial_coeff_parity(unsigned int n, unsigned int k);

/**
 * Computes factorial
 * \param n input value
 * \return n!
 */
unsigned int factorial(unsigned int x, unsigned int result = 1);

/**
 * Computes binomial coefficient
 * \param n n as in (n k)
 * \param k k as in (n k)
 * \return (n k) or 0 if invalid (n<k)
 */
unsigned int binomial_coeff(unsigned int n, unsigned int k);

}
}

#endif // __GF_UTILS_H__
