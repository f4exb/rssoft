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

	 Galois Field GF(q=2^m) class

*/
/*
 *********************************************************************
 *                                                                   *
 *        Galois Field Arithmetic Library (version 0.0.1)            *
 *                                                                   *
 * Class: Galois Field                                               *
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

#ifndef __GFQ_H__
#define __GFQ_H__

#include "GF2_Polynomial.h"
#include <iostream>
#include <vector>
#include <string.h>

namespace rssoft
{
namespace gf
{
typedef int GFq_Symbol;        //!< Symbol or binary-polynomial representation (ex: 5 is X^2+1)
const GFq_Symbol GFERROR = -1; //!< Undefined symbol

/**
 * \brief Galois Field GF(q=2^m) class.
 * Generates and holds lookup tables (LUT) for basic operations.
 * Hosts basic operations and element representation conversions
 */
class GFq
{

public:
	GFq(const int pwr, const GF2_Polynomial& primitive_poly);
	GFq(const GFq& gf);
	~GFq();

	GFq& operator=(const GFq& gf);
	bool operator==(const GFq& gf) const;
	bool operator!=(const GFq& gf) const;

	/**
	 * Alpha based log
	 * \param value Symbol representation of the element
	 * \return The power of alpha corresponding to the element
	 */
	inline unsigned int index(const GFq_Symbol value) const
	{
		return index_of[value];
	}

	/**
	 * Alpha exponentiation
	 * \param value The power of alpha
	 * \return The symbol representation to alpha raised to the given power
	 */
	inline GFq_Symbol alpha(const unsigned int power) const
	{
		return alpha_to[power];
	}

	/**
	 * Get the size of the field
	 * \return The number of elements in the field minus one (i.e. non null elements)
	 */
	inline unsigned int size() const
	{
		return field_size;
	}

	/**
	 * Get m the power of 2 in GF(2^m)
	 * \return GF(2^m) exponent
	 */
	inline unsigned int pwr() const
	{
		return power;
	}

	inline GFq_Symbol add(const GFq_Symbol& a, const GFq_Symbol& b) const
	{
		return (a ^ b);
	}

	inline GFq_Symbol sub(const GFq_Symbol& a, const GFq_Symbol& b) const
	{
		return (a ^ b);
	}

	inline GFq_Symbol mul(const GFq_Symbol& a, const GFq_Symbol& b) const
	{
#if !defined(NO_GFLUT)
		return mul_table[a][b];
#else
		if ((a == 0) || (b == 0))
		{
			return 0;
		}
		else
		{
			return alpha_to[fast_modulus(index_of[a] + index_of[b])];
		}
#endif
	}

	inline GFq_Symbol div(const GFq_Symbol& a, const GFq_Symbol& b) const
	{
#if !defined(NO_GFLUT)
		return div_table[a][b];
#else
		if ((a == 0) || (b == 0))
		{
			return 0;
		}
		else
		{
			return alpha_to[fast_modulus(index_of[a] - index_of[b] + field_size)];
		}
#endif
	}

	inline GFq_Symbol exp(const GFq_Symbol& a, const int& n) const
	{
		if (n == 0)
		{
			return 1;
		}
		else if (a == 0)
        {
            return 0;
        }
        else
        {
        	unsigned int log_a = index_of[a];
        	unsigned int log_a_pwr_n = log_a * n;
        	return alpha_to[log_a_pwr_n % field_size];
        }
	}

/* Just does not work
	inline GFq_Symbol exp(const GFq_Symbol& a, const int& n) const
	{
        if (n == 0)
        {
            return 1;
        }
#if !defined(NO_GFLUT)
		if (n < 0)
		{
			int b = n;
			while (b < 0)
			{
				b += field_size; // b could be negative
			}

			if (b == 0)
			{
				return 1;
			}
			return exp_table[a][b];
		}
		else
		{
			return exp_table[a][n & field_size];
		}
#else
		if (a != 0)
		{
			if (n < 0)
			{
				int b = n;
				while(b < 0)
				{
					b += field_size; // b could be negative
				}
				if (b == 0)
				{
					return 1;
				}
				return alpha_to[fast_modulus(index_of[a] * b)];
			}
			else if (n == 0)
			{
				return 1;
			}
			else
			{
				return alpha_to[fast_modulus(index_of[a] * n)];
			}
		}
		else
		{
			return 0;
		}
#endif
	}
*/

	inline GFq_Symbol inverse(const GFq_Symbol& val) const
	{
#if !defined(NO_GFLUT)
		return mul_inverse[val];
#else
		return alpha_to[fast_modulus(field_size - index_of[val])];
#endif
	}

	friend std::ostream& operator <<(std::ostream& os, const GFq& gf);

private:

	void generate_field();
	GFq_Symbol fast_modulus(GFq_Symbol x) const;
	GFq_Symbol gen_mul(const GFq_Symbol& a, const GFq_Symbol& b) const;
	GFq_Symbol gen_div(const GFq_Symbol& a, const GFq_Symbol& b) const;
	GFq_Symbol gen_exp(const GFq_Symbol& a, const unsigned int& n) const;
	GFq_Symbol gen_inverse(const GFq_Symbol& val) const;

	unsigned int power;                   //!< m the power of 2 as in GF(2^m)
	unsigned int field_size;              //!< Number of non null elements in the field
    const GF2_Polynomial& primitive_poly; //!< Primitive polynomial
	unsigned int prim_poly_hash;          //!< Primitive polynomial hash coded
	GFq_Symbol* alpha_to;                 //!< Exponential or anti-log unary operation LUT
	GFq_Symbol* index_of;                 //!< Log unary operation LUT
	GFq_Symbol* mul_inverse;              //!< Multiplicative inverse unary operation LUT
	GFq_Symbol** mul_table;               //!< Multiplication binary operation LUT
	GFq_Symbol** div_table;               //!< Division binary operation LUT
	GFq_Symbol** exp_table;               //!< Exponent binary operation LUT

};

} // namespace gf
} // namespace rssoft

#endif // __GFQ_H__
