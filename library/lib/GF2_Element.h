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
#ifndef __GF2_ELEMENT_H__
#define __GF2_ELEMENT_H__

#include "GF_Exception.h"
#include <iostream>

namespace rssoft
{
namespace gf
{
typedef unsigned char GF2_Symbol; //!< External representation of symbol in GF(2)

class GF2_Element
{

public:

	GF2_Element(GF2_Symbol v = 0);
	GF2_Element(const GF2_Element& gfe);
	~GF2_Element()
	{
	}

	inline unsigned int uint_value() const
	{
		return (unsigned int) bin_value;
	}

	inline GF2_Element& operator=(const GF2_Element& gfe)
	{
		if (this != &gfe)
		{
			bin_value = gfe.bin_value;
		}

		return *this;
	}

	inline GF2_Element& operator=(const GF2_Symbol& v)
	{
		bin_value = (v ? 1 : 0);
		return *this;
	}

	inline GF2_Element& operator+=(const GF2_Element& gfe)
	{
		bin_value ^= gfe.bin_value;
		return *this;
	}

	inline GF2_Element& operator+=(const GF2_Symbol& v)
	{
		bin_value ^= (v ? 1 : 0);
		return *this;
	}

	inline GF2_Element& operator-=(const GF2_Element& gfe)
	{
		*this += gfe;
		return *this;
	}

	inline GF2_Element& operator-=(const GF2_Symbol& v)
	{
		*this += (v ? 1 : 0);
		return *this;
	}

	inline GF2_Element& operator*=(const GF2_Element& gfe)
	{
		bin_value *= gfe.bin_value;
		return *this;
	}

	inline GF2_Element& operator*=(const GF2_Symbol& v)
	{
		bin_value *= (v ? 1 : 0);
		return *this;
	}

	inline GF2_Element& operator/=(const GF2_Element& gfe)
	{
		if (gfe.bin_value == 0)
		{
			throw GF_Exception("Division by zero");
		}

		return *this;
	}

	inline GF2_Element& operator/=(const GF2_Symbol& v)
	{
		if (!v)
		{
			throw GF_Exception("Division by zero");
		}

		return *this;
	}

	inline GF2_Element& operator^=(const int& n)
	{
		return *this;
	}

	inline bool operator==(const GF2_Element& gfe) const
	{
		return (bin_value == gfe.bin_value);
	}

	inline bool operator==(const GF2_Symbol& v) const
	{
		return (bin_value == (v ? 1 : 0));
	}

	inline bool operator!=(const GF2_Element& gfe) const
	{
		return (bin_value != gfe.bin_value);
	}

	inline bool operator!=(const GF2_Symbol& v) const
	{
		return (bin_value != (v ? 1 : 0));
	}

	inline bool operator<(const GF2_Element& gfe)
	{
		return (bin_value < gfe.bin_value);
	}

	inline bool operator<(const GF2_Symbol& v)
	{
		return (bin_value < (v ? 1 : 0));
	}

	inline bool operator>(const GF2_Element& gfe)
	{
		return (bin_value > gfe.bin_value);
	}

	inline bool operator>(const GF2_Symbol& v)
	{
		return (bin_value > (v ? 1 : 0));
	}

	friend std::ostream& operator <<(std::ostream& os,
			const GF2_Element& gfe);

private:
	unsigned char bin_value; //!< internal representation of the symbol
};

GF2_Element operator +(const GF2_Element& a, const GF2_Element& b);
GF2_Element operator -(const GF2_Element& a, const GF2_Element& b);
GF2_Element operator *(const GF2_Element& a, const GF2_Element& b);
GF2_Element operator *(const GF2_Element& a, const GF2_Symbol& b);
GF2_Element operator *(const GF2_Symbol& a, const GF2_Element& b);
GF2_Element operator /(const GF2_Element& a, const GF2_Element& b);
GF2_Element operator ^(const GF2_Element& a, const int& b);

} // namespacegf
} // namespace rssoft

#endif // __GF2_ELEMENT_H__
