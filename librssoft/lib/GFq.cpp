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

	 Galois Field GF(2^m) class

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

#include "GFq.h"
#include "GF_Exception.h"
#include "GF2_Polynomial.h"

namespace rssoft
{
namespace gf
{

// ================================================================================================
GFq::GFq(const int pwr, const GF2_Polynomial& _primitive_poly) :
		power(pwr), field_size((1 << power) - 1), primitive_poly(_primitive_poly)
{
	if (primitive(primitive_poly, pwr))
	{
			alpha_to = new GFq_Symbol[field_size + 1];
			index_of = new GFq_Symbol[field_size + 1];

		#if !defined(NO_GFLUT)

			mul_table = new GFq_Symbol*[(field_size + 1)];
			div_table = new GFq_Symbol*[(field_size + 1)];
			exp_table = new GFq_Symbol*[(field_size + 1)];
			mul_inverse = new GFq_Symbol[(field_size + 1) * 2];

			for (unsigned int i = 0; i < (field_size + 1); i++)
			{
				mul_table[i] = new GFq_Symbol[(field_size + 1)];
				div_table[i] = new GFq_Symbol[(field_size + 1)];
				exp_table[i] = new GFq_Symbol[(field_size + 1)];
			}

		#else

			mul_table = new GFq_Symbol *[1];
			div_table = new GFq_Symbol *[1];
			exp_table = new GFq_Symbol *[1];
			mul_inverse = new GFq_Symbol [1];

		#endif

			prim_poly_hash = 0xAAAAAAAA;
			unsigned int prim_poly_mono_power = 1;

			for (unsigned int i = 0; i <= power; i++)
			{
				if (i < power)
				{
					prim_poly_hash += ((i & 1) == 0) ? (  (prim_poly_hash <<  7) ^ primitive_poly[i].uint_value() ^ (prim_poly_hash >> 3)) :
									  (~((prim_poly_hash << 11) ^ primitive_poly[i].uint_value() ^ (prim_poly_hash >> 5)));
				}

				prim_poly_mono_power <<= 1;
			}

			generate_field();
	}
	else
	{
		throw GF_Exception("Non primitive polynomial used to create GF(q) field");
	}
}


// ================================================================================================
GFq::GFq(const GFq& gf) :
		primitive_poly(gf.primitive_poly)
{
	power = gf.power;
	field_size = gf.field_size;
	prim_poly_hash = gf.prim_poly_hash;
	alpha_to = new GFq_Symbol[field_size + 1];
	index_of = new GFq_Symbol[field_size + 1];

	memcpy(alpha_to, gf.alpha_to, (field_size + 1) * sizeof(GFq_Symbol));
	memcpy(index_of, gf.index_of, (field_size + 1) * sizeof(GFq_Symbol));

#if !defined(NO_GFLUT)

	mul_table = new GFq_Symbol*[(field_size + 1)];
	div_table = new GFq_Symbol*[(field_size + 1)];
	exp_table = new GFq_Symbol*[(field_size + 1)];
	mul_inverse = new GFq_Symbol[(field_size + 1) * 2];

	for (unsigned int i = 0; i < (field_size + 1); i++)
	{
		mul_table[i] = new GFq_Symbol[(field_size + 1)];
		div_table[i] = new GFq_Symbol[(field_size + 1)];
		exp_table[i] = new GFq_Symbol[(field_size + 1)];
	}

	memcpy(mul_inverse, gf.mul_inverse, (field_size + 1) * sizeof(GFq_Symbol) * 2);

	memcpy(mul_table, gf.mul_table, (field_size + 1) * sizeof(GFq_Symbol*));
	memcpy(div_table, gf.div_table, (field_size + 1) * sizeof(GFq_Symbol*));
	memcpy(exp_table, gf.exp_table, (field_size + 1) * sizeof(GFq_Symbol*));

	for (unsigned int i = 0; i < (field_size + 1); i++)
	{
		memcpy(mul_table[i], gf.mul_table[i], (field_size + 1) * sizeof(GFq_Symbol));
		memcpy(div_table[i], gf.div_table[i], (field_size + 1) * sizeof(GFq_Symbol));
		memcpy(exp_table[i], gf.exp_table[i], (field_size + 1) * sizeof(GFq_Symbol));
	}

#endif
}


// ================================================================================================
GFq::~GFq()
{
	if (alpha_to != NULL)
	{
		delete[] alpha_to;
	}
	if (index_of != NULL)
	{
		delete[] index_of;
	}

#if !defined(NO_GFLUT)

	for (unsigned int i = 0; i < (field_size + 1); i++)
	{
		if (mul_table[i] != NULL)
		{
			delete[] mul_table[i];
		}
		if (div_table[i] != NULL)
		{
			delete[] div_table[i];
		}
		if (exp_table[i] != NULL)
		{
			delete[] exp_table[i];
		}
	}

	if (mul_table != NULL)
	{
		delete[] mul_table;
	}

	if (div_table != NULL)
	{
		delete[] div_table;
	}

	if (exp_table != NULL)
	{
		delete[] exp_table;
	}

	if (mul_inverse != NULL)
	{
		delete[] mul_inverse;
	}

#endif
}


// ================================================================================================
bool GFq::operator==(const GFq& gf) const
{
	return ((this->power == gf.power) && (this->prim_poly_hash == gf.prim_poly_hash));
}


// ================================================================================================
bool GFq::operator!=(const GFq& gf) const
{
	return ((this->power != gf.power) || (this->prim_poly_hash != gf.prim_poly_hash));
}


// ================================================================================================
GFq& GFq::operator=(const GFq& gf)
{
	if (this == &gf)
	{
		return *this;
	}

	if (alpha_to != NULL)
	{
		delete[] alpha_to;
	}

	if (index_of != NULL)
	{
		delete[] index_of;
	}

	power = gf.power;
	field_size = gf.field_size;
	prim_poly_hash = gf.prim_poly_hash;

	memcpy(alpha_to, gf.alpha_to, (field_size + 1) * sizeof(GFq_Symbol));
	memcpy(index_of, gf.index_of, (field_size + 1) * sizeof(GFq_Symbol));

#if !defined(NO_GFLUT)

	if (mul_table != NULL)
	{
		delete[] mul_table;
	}

	if (div_table != NULL)
	{
		delete[] div_table;
	}

	if (exp_table != NULL)
	{
		delete[] exp_table;
	}

	if (mul_inverse != NULL)
	{
		delete[] mul_inverse;
	}

	mul_table = new GFq_Symbol*[(field_size + 1)];
	div_table = new GFq_Symbol*[(field_size + 1)];
	exp_table = new GFq_Symbol*[(field_size + 1)];
	mul_inverse = new GFq_Symbol[(field_size + 1) * 2];

	for (unsigned int i = 0; i < (field_size + 1); i++)
	{
		mul_table[i] = new GFq_Symbol[(field_size + 1)];
		div_table[i] = new GFq_Symbol[(field_size + 1)];
		exp_table[i] = new GFq_Symbol[(field_size + 1)];
	}

	memcpy(mul_inverse, gf.mul_inverse, (field_size + 1) * sizeof(GFq_Symbol) * 2);

	memcpy(mul_table, gf.mul_table, (field_size + 1) * sizeof(GFq_Symbol*));
	memcpy(div_table, gf.div_table, (field_size + 1) * sizeof(GFq_Symbol*));
	memcpy(exp_table, gf.exp_table, (field_size + 1) * sizeof(GFq_Symbol*));

	for (unsigned int i = 0; i < (field_size + 1); i++)
	{
		memcpy(mul_table[i], gf.mul_table[i], (field_size + 1) * sizeof(GFq_Symbol));
		memcpy(div_table[i], gf.div_table[i], (field_size + 1) * sizeof(GFq_Symbol));
		memcpy(exp_table[i], gf.exp_table[i], (field_size + 1) * sizeof(GFq_Symbol));
	}

#endif

	return *this;
}


// ================================================================================================
void GFq::generate_field()
{
	/*
	 Note: It is assumed that the degree of the primitive
	 polynomial will be equivelent to the m value as
	 in GF(2^m)
	 */

	/*
	 need to update using stanford method for prim-poly generation.
	 */
	int mask = 1;

	alpha_to[power] = 0;

	for (unsigned int i = 0; i < power; i++)
	{
		alpha_to[i] = mask;
		index_of[alpha_to[i]] = i;

		if (primitive_poly[i] != 0)
		{
			alpha_to[power] ^= mask;
		}

		mask <<= 1;
	}

	index_of[alpha_to[power]] = power;

	mask >>= 1;

	for (unsigned int i = power + 1; i < field_size; i++)
	{
		if (alpha_to[i - 1] >= mask)
		{
			alpha_to[i] = alpha_to[power] ^ ((alpha_to[i - 1] ^ mask) << 1);
		}
		else
		{
			alpha_to[i] = alpha_to[i - 1] << 1;
		}

		index_of[alpha_to[i]] = i;
	}

	index_of[0] = GFERROR;
	alpha_to[field_size] = 1;

#if !defined(NO_GFLUT)

	for (unsigned int i = 0; i < field_size + 1; i++)
	{
		for (unsigned int j = 0; j < field_size + 1; j++)
		{
			mul_table[i][j] = gen_mul(i, j);
			div_table[i][j] = gen_div(i, j);
			exp_table[i][j] = gen_exp(i, j);
		}
	}

	for (unsigned int i = 0; i < (field_size + 1); i++)
	{
		mul_inverse[i] = gen_inverse(i);
		mul_inverse[i + (field_size + 1)] = mul_inverse[i];
	}

#endif
}


// ================================================================================================
GFq_Symbol GFq::fast_modulus(GFq_Symbol x) const
{
	while (x >= (int) field_size)
	{
		x -= (int) field_size;
		x = (x >> power) + (x & (int) field_size);
	}

	return x;
}


// ================================================================================================
GFq_Symbol GFq::gen_mul(const GFq_Symbol& a, const GFq_Symbol& b) const
{
	if ((a == 0) || (b == 0))
	{
		return 0;
	}
	else
	{
		return alpha_to[fast_modulus(index_of[a] + index_of[b])];
	}
}


// ================================================================================================
GFq_Symbol GFq::gen_div(const GFq_Symbol& a, const GFq_Symbol& b) const
{
	if ((a == 0) || (b == 0))
	{
		return 0;
	}
	else
	{
		return alpha_to[fast_modulus(index_of[a] - index_of[b] + field_size)];
	}
}


// ================================================================================================
GFq_Symbol GFq::gen_exp(const GFq_Symbol& a, const unsigned int& n) const
{
	if (a != 0)
	{
		if (n == 0)
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
}


// ================================================================================================
GFq_Symbol GFq::gen_inverse(const GFq_Symbol& val) const
{
	return alpha_to[fast_modulus(field_size - index_of[val])];
}


// ================================================================================================
std::ostream& operator << (std::ostream& os, const GFq& gf)
{
    os << "GF(2^" << gf.pwr() << ")" << std::endl;

    os << "P = " << gf.primitive_poly << std::endl;
    os << "i\ta^i\tlog_a(i)" << std::endl;

    for(unsigned int i = 0; i < gf.field_size + 1; i++)
    {
        os << i << "\t" << gf.alpha_to[i] << "\t" << gf.index_of[i] << std::endl;
    }

    return os;
}

} // namespace gf
} // namespace rssoft
