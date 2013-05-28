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

 Univariate polynomials with coefficients in GF(2) class hence
 members of GF(2)[X]

 */

#include "GF2_Polynomial.h"
#include "GF_Exception.h"
#include <algorithm>

namespace rssoft
{
namespace gf
{

// ================================================================================================
GF2_Polynomial::GF2_Polynomial()
{
	poly.clear();
}

// ================================================================================================
GF2_Polynomial::GF2_Polynomial(unsigned int size, GF2_Element* gfe)
{
	if (gfe != NULL)
	{
		poly.assign(gfe, gfe + size);
	}
	else
	{
		poly.assign(size, 0);
	}
}

// ================================================================================================
GF2_Polynomial::GF2_Polynomial(const GF2_Polynomial& polynomial)
{
	poly = polynomial.poly;
}

// ================================================================================================
GF2_Polynomial::GF2_Polynomial(const GF2_Element& gfe)
{
	poly.clear();
	poly.push_back(gfe);
}

// ================================================================================================
GF2_Polynomial::GF2_Polynomial(unsigned int n)
{
	poly.assign(n+1, 0);
	poly.back() = 1;
}

// ================================================================================================
bool GF2_Polynomial::valid() const
{
	return (poly.size() > 0);
}

// ================================================================================================
bool GF2_Polynomial::null() const
{
	return (poly.size() == 0) || ((poly.size() == 1) && (poly[0] == 0));
}

// ================================================================================================
unsigned int GF2_Polynomial::deg() const
{
	return static_cast<unsigned int>(poly.size() - 1);
}

// ================================================================================================
const std::vector<GF2_Element>& GF2_Polynomial::get_poly() const
{
	return poly;
}

// ================================================================================================
std::vector<GF2_Element>& GF2_Polynomial::get_poly_for_update()
{
	return poly;
}

// ================================================================================================
void GF2_Polynomial::set_degree(unsigned int x)
{
	poly.resize(x - 1, 0);
}

// ================================================================================================
GF2_Polynomial& GF2_Polynomial::operator=(const GF2_Polynomial& polynomial)
{
	if (this != &polynomial)
	{
		poly = polynomial.poly;
	}
	return *this;
}

// ================================================================================================
GF2_Polynomial& GF2_Polynomial::operator=(const GF2_Element& gfe)
{
	poly.clear();
	poly.push_back(gfe);
	return *this;
}

// ================================================================================================
GF2_Polynomial& GF2_Polynomial::operator+=(const GF2_Polynomial& polynomial)
{
	if (poly.size() < polynomial.poly.size())
	{
		unsigned int j = 0;
		for (unsigned int i = 0; i < poly.size(); i++)
		{
			poly[i] += polynomial.poly[j++];
		}

		for (; j < polynomial.poly.size(); j++)
		{
			poly.push_back(polynomial.poly[j]);
		}
	}
	else
	{
		unsigned int i = 0;
		for (unsigned int j = 0; j < polynomial.poly.size(); j++)
		{
			poly[i++] += polynomial.poly[j];
		}
	}

	simplify(*this);

	return *this;
}

// ================================================================================================
GF2_Polynomial& GF2_Polynomial::operator+=(const GF2_Element& gfe)
{
	poly[0] += gfe;
	return *this;
}

// ================================================================================================
GF2_Polynomial& GF2_Polynomial::operator-=(const GF2_Polynomial& polynomial)
{
	return (*this += polynomial);
}

// ================================================================================================
GF2_Polynomial& GF2_Polynomial::operator-=(const GF2_Element& gfe)
{
	poly[0] -= gfe;
	return *this;
}

// ================================================================================================
GF2_Polynomial& GF2_Polynomial::operator*=(const GF2_Polynomial& polynomial)
{
	GF2_Polynomial product(deg() + polynomial.deg() + 1, 0); // create polynomial with all zero coefficients

	for (unsigned int i = 0; i < poly.size(); i++)
	{
		for (unsigned int j = 0; j < polynomial.poly.size(); j++)
		{
			product.poly[i + j] += poly[i] * polynomial.poly[j];
		}
	}

	simplify(product);
	poly = product.poly;

	return *this;
}

// ================================================================================================
GF2_Polynomial& GF2_Polynomial::operator*=(const GF2_Element& gfe)
{
	if (gfe == 0)
	{
		poly.clear();
		poly.push_back(0);
	}

	return *this;
}

// ================================================================================================
GF2_Polynomial& GF2_Polynomial::operator/=(const GF2_Polynomial& divisor)
{
	const std::pair<GF2_Polynomial, GF2_Polynomial>& divres = div(*this, divisor);
	poly = divres.first.get_poly();
	return *this;
}

// ================================================================================================
GF2_Polynomial& GF2_Polynomial::operator/=(const GF2_Element& gfe)
{
	if (gfe == 0)
	{
		throw GF_Exception("Division by zero");
	}

	return *this;
}

// ================================================================================================
GF2_Polynomial& GF2_Polynomial::operator%=(const GF2_Polynomial& divisor)
{
	const std::pair<GF2_Polynomial, GF2_Polynomial>& divres = div(*this, divisor);
	poly = divres.second.get_poly();
	return *this;
}

// ================================================================================================
GF2_Polynomial& GF2_Polynomial::operator%=(const unsigned int& power)
{
	if (poly.size() >= power)
	{
		for (size_t i = poly.size(); i > power; i--)
		{
			poly.pop_back();
		}
	}

	simplify(*this);
	return *this;
}

// ================================================================================================
GF2_Polynomial& GF2_Polynomial::operator^=(const int& n)
{
	GF2_Polynomial result = *this;

	for (int i = 0; i < n; i++)
	{
		result *= *this;
	}

	*this = result;

	return *this;
}

// ================================================================================================
GF2_Polynomial& GF2_Polynomial::operator<<=(const unsigned int& n)
{
	if (poly.size() > 0)
	{
		std::size_t initial_size = poly.size();
		poly.resize(poly.size() + n, 1);
		std::copy(poly.rend() - initial_size, poly.rend(), poly.rbegin());
		std::fill(poly.begin(), poly.begin() + n, 0);
	}

	return *this;
}

// ================================================================================================
GF2_Polynomial& GF2_Polynomial::operator>>=(const unsigned int& n)
{
	if (n <= poly.size())
	{
		std::copy(poly.begin() + n, poly.end(), poly.begin());
		poly.erase(poly.end() - n, poly.end());
	}
	else
	{
		poly.clear();
	}
	return *this;
}

// ================================================================================================
GF2_Element& GF2_Polynomial::operator[](unsigned int term)
{
	assert(term < poly.size());
	return poly[term];
}

// ================================================================================================
const GF2_Element& GF2_Polynomial::operator[](unsigned int term) const
{
	assert(term < poly.size());
	return poly[term];
}

// ================================================================================================
GF2_Element GF2_Polynomial::operator()(GF2_Element value)
{
	GF2_Element result = 0;

	if (poly.size() > 0)
	{
		result = poly[poly.size() - 1];
		for (std::size_t i = poly.size() - 2; ((int) i) >= 0; i--)
		{
			result = poly[i] + (result * value);
		}
		return result;
	}
	else
	{
		throw GF_Exception("Cannot evaluate invalid polynomial");
	}
}

// ================================================================================================
const GF2_Element GF2_Polynomial::operator()(GF2_Element value) const
{
	GF2_Element result = 0;

	if (poly.size() > 0)
	{
		result = poly[poly.size() - 1];
		for (std::size_t i = poly.size() - 2; ((int) i) >= 0; i--)
		{
			result = poly[i] + (result * value);
		}
		return result;
	}
	else
	{
		throw GF_Exception("Cannot evaluate invalid polynomial");
	}
}

// ================================================================================================
bool GF2_Polynomial::operator==(const GF2_Polynomial& polynomial) const
{
	if (poly.size() != polynomial.poly.size())
		return false;
	else
	{
		for (unsigned int i = 0; i < poly.size(); i++)
		{
			if (poly[i] != polynomial.poly[i])
				return false;
		}
		return true;
	}
}

// ================================================================================================
bool GF2_Polynomial::operator!=(const GF2_Polynomial& polynomial) const
{
	return !(*this == polynomial);
}

// ================================================================================================
void simplify(GF2_Polynomial& polynomial)
{
	if (polynomial.get_poly().size() > 0)
	{
		for (size_t i = polynomial.get_poly().size() - 1; i > 0; i--)
		{
			if (polynomial.get_poly()[i] == 0)
			{
				polynomial.get_poly_for_update().pop_back();
			}
			else
			{
				break;
			}
		}
	}
}

// ================================================================================================
GF2_Polynomial operator+(const GF2_Polynomial& a, const GF2_Polynomial& b)
{
	GF2_Polynomial result = a;
	result += b;
	return result;
}

// ================================================================================================
GF2_Polynomial operator +(const GF2_Polynomial& a, const GF2_Element& b)
{
	GF2_Polynomial result = a;
	result += b;
	return result;
}

// ================================================================================================
GF2_Polynomial operator +(const GF2_Element& a, const GF2_Polynomial& b)
{
	GF2_Polynomial result = b;
	result += a;
	return result;
}

// ================================================================================================
GF2_Polynomial operator -(const GF2_Polynomial& a, const GF2_Polynomial& b)
{
	GF2_Polynomial result = a;
	result -= b;
	return result;
}

// ================================================================================================
GF2_Polynomial operator -(const GF2_Polynomial& a, const GF2_Element& b)
{
	GF2_Polynomial result = a;
	result -= b;
	return result;
}

// ================================================================================================
GF2_Polynomial operator -(const GF2_Element& a, const GF2_Polynomial& b)
{
	GF2_Polynomial result = b;
	result -= a;
	return result;
}

// ================================================================================================
GF2_Polynomial operator *(const GF2_Polynomial& a, const GF2_Polynomial& b)
{
	GF2_Polynomial result = a;
	result *= b;
	return result;
}

// ================================================================================================
GF2_Polynomial operator *(const GF2_Element& a, const GF2_Polynomial& b)
{
	GF2_Polynomial result = b;
	result *= a;
	return result;
}

// ================================================================================================
GF2_Polynomial operator *(const GF2_Polynomial& a, const GF2_Element& b)
{
	GF2_Polynomial result = a;
	result *= b;
	return result;
}

// ================================================================================================
GF2_Polynomial operator /(const GF2_Polynomial& a, const GF2_Polynomial& b)
{
	GF2_Polynomial result = a;
	result /= b;
	return result;
}

// ================================================================================================
GF2_Polynomial operator /(const GF2_Polynomial& a, const GF2_Element& b)
{
	GF2_Polynomial result = a;
	result /= b;
	return result;
}

// ================================================================================================
GF2_Polynomial operator %(const GF2_Polynomial& a, const GF2_Polynomial& b)
{
	GF2_Polynomial result = a;
	result %= b;
	return result;
}

// ================================================================================================
GF2_Polynomial operator %(const GF2_Polynomial& a, const unsigned int& power)
{
	GF2_Polynomial result = a;
	result %= power;
	return result;
}

// ================================================================================================
GF2_Polynomial operator ^(const GF2_Polynomial& a, const int& n)
{
	GF2_Polynomial result = a;
	result ^= n;
	return result;
}

// ================================================================================================
GF2_Polynomial operator<<(const GF2_Polynomial& a, const unsigned int& n)
{
	GF2_Polynomial result = a;
	result <<= n;
	return result;
}

// ================================================================================================
GF2_Polynomial operator>>(const GF2_Polynomial& a, const unsigned int& n)
{
	GF2_Polynomial result = a;
	result >>= n;
	return result;
}

// ================================================================================================
GF2_Polynomial gcd(const GF2_Polynomial& a, const GF2_Polynomial& b)
{
	if (a.null() && b.null())
	{
		throw GF_Exception("GCD with both zero operand polynomials");
	}

	if (a.null())
	{
		return b;
	}

	if (b.null())
	{
		return a;
	}

	GF2_Polynomial r = ((a.deg() < b.deg()) ? a : b);
	GF2_Polynomial x = ((a.deg() < b.deg()) ? b : a);
	GF2_Polynomial t = ((a.deg() < b.deg()) ? b : a);

	while (!r.null())
	{
		t = r;
		r = x % t;
		x = t;
	}

	return x;
}

// ================================================================================================
std::pair<GF2_Polynomial, GF2_Polynomial> div(const GF2_Polynomial& dividend,
		const GF2_Polynomial& divisor)
{
	if (divisor.deg() < 0)
	{
		throw GF_Exception("GF Poly Division invalid operands");
	}
	else
	{
		if (divisor.deg() == 0)
		{
			GF2_Polynomial quotient(dividend);
			quotient /= divisor[0];
			return std::make_pair(quotient, GF2_Polynomial((GF2_Element)0));
		}
		else if (dividend.deg() < divisor.deg())
		{
			return std::make_pair(dividend, GF2_Polynomial((GF2_Element)0));
		}
		else
		{
			GF2_Polynomial remainder(dividend);
			GF2_Polynomial quotient(dividend.deg() - divisor.deg() + 1, 0); // creates a polynomial with all zero coefficients

			while (remainder.valid() && (remainder.deg() >= divisor.deg()))
			{
				quotient[remainder.deg() - divisor.deg()] = remainder[remainder.deg()] / divisor[divisor.deg()];

				int r,d;

				for (r=remainder.deg(), d=divisor.deg(); d >=0; r--, d--)
				{
					remainder[r] -= quotient[remainder.deg() - divisor.deg()]*divisor[d];
				}

				simplify(remainder);
			}

			simplify(quotient);

			return std::make_pair(quotient, remainder);
		}
	}
}

// ================================================================================================
bool irreducible(const GF2_Polynomial& f)
{
	int df = f.deg();
	
	// trivial answers
	if (df <= 0) 
	{
		return false;
	}
	else if (df == 1)
	{
		return true;
	}
	// real stuff...
	else
	{
		rssoft::gf::GF2_Element Te[3] = {
			rssoft::gf::GF2_Element(0),
			rssoft::gf::GF2_Element(1),
			rssoft::gf::GF2_Element(1),
		};
		size_t Te_sz = sizeof(Te)/sizeof(rssoft::gf::GF2_Element);
		rssoft::gf::GF2_Polynomial T(Te_sz,Te); // T(x) = x^2 + x
		
		rssoft::gf::GF2_Polynomial One(rssoft::gf::GF2_Element(1)); // One(x)=1
		
		return (gcd(f,T) == One); // or gcd(f,(X2-X)%f) but for deg >= 2 this does not matter
	}
}

// ================================================================================================
unsigned int coeff_parity(const GF2_Polynomial& a)
{
	const std::vector<GF2_Element>& poly = a.get_poly();
	std::vector<GF2_Element>::const_iterator it = poly.begin();
	unsigned int count=0;
	
	for (; it != poly.end(); ++it)
	{
		if (*it != 0)
		{
			count++;
		}
	}
	
	return count;
}

// ================================================================================================
bool primitive(const GF2_Polynomial& a, unsigned int m)
{
	unsigned int k = 1;
	k <<= m;
	k--; // 2^m-1
	
	rssoft::gf::GF2_Polynomial One(rssoft::gf::GF2_Element(1)); // One(x)=1
	rssoft::gf::GF2_Polynomial Xk(k); // Xk(x) = x^k
	
	if (irreducible(a) && (a.deg() == m) && ((coeff_parity(a)%2)==1))
	{
		if (((Xk+One)%a).null())
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}

// ================================================================================================
std::ostream& operator <<(std::ostream& os, const GF2_Polynomial& polynomial)
{
	if (polynomial.deg() >= 0)
	{
		bool is_null = true;
		bool first_coeff = true;

		for (unsigned int i = 0; i < polynomial.poly.size(); i++)
		{
			GF2_Element coeff = polynomial.poly[i];

			if (coeff != 0)
			{
				is_null = false;
				os << ((!first_coeff) ? "+ " : "");

				if (i == 0)
				{
					os << "1 ";
				}
				else if (i==1)
				{
					os << "x ";
				}
				else
				{
					os << "x^" << i << " ";
				}

				first_coeff = false;
			}
		}

		if (is_null)
		{
			os << "0";
		}
	}

	return os;
}

} // namespace gf
} // namespace rssoft
