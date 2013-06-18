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
 * Edouard Griffiths, 2013: modified to secure constructor           *
 *                                                                   *
 * Copyright Notice:                                                 *
 * Free use of this library is permitted under the guidelines and    *
 * in accordance with the most current version of the Common Public  *
 * License.                                                          *
 * http://www.opensource.org/licenses/cpl.php                        *
 *                                                                   *
 *********************************************************************
 */

#include "GFq_Element.h"

namespace rssoft
{
namespace gf
{

	GFq_Element::GFq_Element(const GFq& _gf, GFq_Symbol v) :
	gf(_gf)
	{
		poly_value = v;
	}

	GFq_Element::GFq_Element(const GFq_Element& gfe) :
	gf(gfe.gf)
	{
		poly_value = gfe.poly_value;
	}

	std::ostream& operator << (std::ostream& os, const GFq_Element& gfe)
	{
		//os << gfe.poly_value;
		if (gfe.poly_value == 0)
		{
			os << "0";
		}
		else if (gfe.poly_value == 1)
		{
			os << "1";
		}
		else
		{
			os << "a^" << gfe.field().index(gfe.poly_value);
		}
		return os;
	}
    
    GFq_Symbol gfq_element_to_symbol(const GFq_Element& gfe)
    {
        return gfe.poly();
    }

	GFq_Element operator+(const GFq_Element& a, const GFq_Element& b)
	{
		GFq_Element result = a;
		result += b;
		return result;
	}

	GFq_Element operator-(const GFq_Element& a, const GFq_Element& b)
	{
		GFq_Element result = a;
		result -= b;
		return result;
	}

	GFq_Element operator*(const GFq_Element& a, const GFq_Element& b)
	{
		GFq_Element result = a;
		result *= b;
		return result;
	}

	GFq_Element operator*(const GFq_Element& a, const GFq_Symbol& b)
	{
		GFq_Element result = a;
		result *= b;
		return result;
	}

	GFq_Element operator*(const GFq_Symbol& a, const GFq_Element& b)
	{
		GFq_Element result = b;
		result *= a;
		return result;
	}

	GFq_Element operator/(const GFq_Element& a, const GFq_Element& b)
	{
		GFq_Element result = a;
		result /= b;
		return result;
	}

	GFq_Element operator^(const GFq_Element& a, const int& b)
	{
		GFq_Element result = a;
		result ^= b;
		return result;
	}

} // namespace gf
} // namespace rssoft

