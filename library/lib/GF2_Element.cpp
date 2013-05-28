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

 Element in Galois Field GF(2) class

 */

#include "GF2_Element.h"

namespace rssoft
{
namespace gf
{

GF2_Element::GF2_Element(GF2_Symbol v)
{
	bin_value = (v ? 1 : 0);
}

GF2_Element::GF2_Element(const GF2_Element& gfe)
{
	bin_value = gfe.bin_value;
}

std::ostream& operator <<(std::ostream& os, const GF2_Element& gfe)
{
	os << (gfe.bin_value ? 1 : 0);
	return os;
}

GF2_Element operator+(const GF2_Element& a, const GF2_Element& b)
{
	GF2_Element result = a;
	result += b;
	return result;
}

GF2_Element operator-(const GF2_Element& a, const GF2_Element& b)
{
	GF2_Element result = a;
	result -= b;
	return result;
}

GF2_Element operator*(const GF2_Element& a, const GF2_Element& b)
{
	GF2_Element result = a;
	result *= b;
	return result;
}

GF2_Element operator*(const GF2_Element& a, const GF2_Symbol& b)
{
	GF2_Element result = a;
	result *= b;
	return result;
}

GF2_Element operator*(const GF2_Symbol& a, const GF2_Element& b)
{
	GF2_Element result = b;
	result *= a;
	return result;
}

GF2_Element operator/(const GF2_Element& a, const GF2_Element& b)
{
	GF2_Element result = a;
	result /= b;
	return result;
}

GF2_Element operator^(const GF2_Element& a, const int& b)
{
	GF2_Element result = a;
	result ^= b;
	return result;
}

} // namespace gf
} // namespace rssoft

