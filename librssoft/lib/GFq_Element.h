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

	 Element in Galois Field GF(q=2^m) class

*/
/*
 *********************************************************************
 *                                                                   *
 *        Galois Field Arithmetic Library (version 0.0.1)            *
 *                                                                   *
 * Class: Galois Field Element                                       *
 * Version: 0.0.1                                                    *
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

#ifndef __GFQ_ELEMENT_H__
#define __GFQ_ELEMENT_H__

#include <iostream>
#include <vector>
#include "GFq.h"

namespace rssoft
{
namespace gf
{

class GFq_Element
{

public:

	GFq_Element(const GFq& _gf, GFq_Symbol v = 0);
	GFq_Element(const GFq_Element& gfe);
	~GFq_Element()
	{
	}

	inline GFq_Element& operator=(const GFq_Element& gfe)
	{
		if (this == &gfe)
			return *this;

		const_cast<GFq&>(gf) = gfe.gf;
		poly_value = gfe.poly_value;

		return *this;
	}

	inline GFq_Element& operator=(const GFq_Symbol& v)
	{
		poly_value = v & gf.size();
		return *this;
	}

	inline GFq_Element& operator+=(const GFq_Element& gfe)
	{
		poly_value ^= gfe.poly_value;
		return *this;
	}

	inline GFq_Element& operator+=(const GFq_Symbol& v)
	{
		poly_value ^= v;
		return *this;
	}

	inline GFq_Element& operator-=(const GFq_Element& gfe)
	{
		*this += gfe;
		return *this;
	}

	inline GFq_Element& operator-=(const GFq_Symbol& v)
	{
		*this += v;
		return *this;
	}

	inline GFq_Element& operator*=(const GFq_Element& gfe)
	{
		poly_value = gf.mul(poly_value, gfe.poly_value);
		return *this;
	}

	inline GFq_Element& operator*=(const GFq_Symbol& v)
	{
		poly_value = gf.mul(poly_value, v);
		return *this;
	}

	inline GFq_Element& operator/=(const GFq_Element& gfe)
	{
		poly_value = gf.div(poly_value, gfe.poly_value);
		return *this;
	}

	inline GFq_Element& operator/=(const GFq_Symbol& v)
	{
		poly_value = gf.div(poly_value, v);
		return *this;
	}

	inline GFq_Element& operator^=(const int& n)
	{
		poly_value = gf.exp(poly_value, n);
		return *this;
	}

	inline bool operator==(const GFq_Element& gfe) const
	{
		return ((gf == gfe.gf) && (poly_value == gfe.poly_value));
	}

	inline bool operator==(const GFq_Symbol& v) const
	{
		return (poly_value == v);
	}

	inline bool operator!=(const GFq_Element& gfe) const
	{
		return ((gf != gfe.gf) || (poly_value != gfe.poly_value));
	}

	inline bool operator!=(const GFq_Symbol& v) const
	{
		return (poly_value != v);
	}

	inline bool operator<(const GFq_Element& gfe) const
	{
		return (poly_value < gfe.poly_value);
	}

	inline bool operator<(const GFq_Symbol& v) const
	{
		return (poly_value < v);
	}

	inline bool operator>(const GFq_Element& gfe) const
	{
		return (poly_value > gfe.poly_value);
	}

	inline bool operator>(const GFq_Symbol& v) const
	{
		return (poly_value > v);
	}

	inline GFq_Symbol index() const
	{
		return gf.index(poly_value);
	}

	inline GFq_Symbol poly() const
	{
		return poly_value;
	}

	inline const GFq& field() const
	{
		return gf;
	}

	inline GFq_Symbol inverse() const
	{
		return gf.inverse(poly_value);
	}

	inline bool is_zero() const
	{
		return poly_value == 0;
	}

	inline bool is_one() const
	{
		return poly_value == 1;
	}
    
	friend std::ostream& operator <<(std::ostream& os,
			const GFq_Element& gfe);
            
private:

	const GFq& gf;
	GFq_Symbol poly_value;
};

GFq_Element operator +(const GFq_Element& a, const GFq_Element& b);
GFq_Element operator -(const GFq_Element& a, const GFq_Element& b);
GFq_Element operator *(const GFq_Element& a, const GFq_Element& b);
GFq_Element operator *(const GFq_Element& a, const GFq_Symbol& b);
GFq_Element operator *(const GFq_Symbol& a, const GFq_Element& b);
GFq_Element operator /(const GFq_Element& a, const GFq_Element& b);
GFq_Element operator ^(const GFq_Element& a, const int& b);

GFq_Symbol gfq_element_to_symbol(const GFq_Element& gfe);

} // namespace gf
} // namespace rssoft

#endif // __GFQ_ELEMENT_H__
