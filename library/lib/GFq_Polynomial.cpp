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
 Extensively modified augmented with a few methods and included in
 RSSoft in gf sub-namespace.

 Univariate polynomials with coefficients in GF(2^m) class hence
 members of GF(2^m)[X]

 */
/*
 *********************************************************************
 *                                                                   *
 *        Galois Field Arithmetic Library (version 0.0.1)            *
 *                                                                   *
 * Class: Galois Field Polynomial                                    *
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

#include "GFq_Polynomial.h"
#include "GF_Exception.h"
#include <algorithm>
#include <numeric>

namespace rssoft
{
namespace gf
{

// ================================================================================================
GFq_Polynomial::GFq_Polynomial(const GFq& _gf) :
		gf(_gf), alpha_format(false)
{
	poly.clear();
}

// ================================================================================================
GFq_Polynomial::GFq_Polynomial(const GFq& _gf, const unsigned int size, GFq_Element* gfe) :
		gf(_gf), alpha_format(false)
{
	if (gfe != NULL)
	{
		poly.assign(gfe, gfe + size);
	}
	else
	{
		poly.assign(size, GFq_Element(gf, 0));
	}
}

// ================================================================================================
GFq_Polynomial::GFq_Polynomial(const GFq& _gf, const std::vector<GFq_Element>& gfe) :
		gf(_gf),
		alpha_format(false),
		poly(gfe.begin(), gfe.end())
{
}

// ================================================================================================
GFq_Polynomial::GFq_Polynomial(const GFq_Polynomial& polynomial) :
		gf(polynomial.gf)
{
	poly = polynomial.poly;
	alpha_format = polynomial.alpha_format;
}

// ================================================================================================
GFq_Polynomial::GFq_Polynomial(const GFq_Element& gfe) :
		gf(gfe.field()), alpha_format(false)
{
	poly.clear();
	poly.push_back(gfe);
}

// ================================================================================================
void GFq_Polynomial::init(const std::vector<GFq_Element>& _poly)
{
	poly = _poly;
}


// ================================================================================================
bool GFq_Polynomial::is_valid() const
{
	return (poly.size() > 0);
}

// ================================================================================================
bool GFq_Polynomial::is_zero() const
{
	return (poly.size() == 0) || ((poly.size() == 1) && (poly[0] == 0));
}

// ================================================================================================
bool GFq_Polynomial::is_one() const
{
	return ((poly.size() == 1) && (poly[0] == 1));
}

// ================================================================================================
unsigned int GFq_Polynomial::deg() const
{
	return static_cast<unsigned int>(poly.size() - 1);
}

// ================================================================================================
const GFq& GFq_Polynomial::field() const
{
	return gf;
}

// ================================================================================================
const std::vector<GFq_Element>& GFq_Polynomial::get_poly() const
{
	return poly;
}

// ================================================================================================
std::vector<GFq_Element>& GFq_Polynomial::get_poly_for_update()
{
	return poly;
}

// ================================================================================================
void GFq_Polynomial::set_degree(const unsigned int& x)
{
	poly.resize(x - 1, GFq_Element(gf, 0));
}

// ================================================================================================
GFq_Polynomial& GFq_Polynomial::operator=(const GFq_Polynomial& polynomial)
{
	if (this == &polynomial)
	{
		return *this;
	}

	const_cast<GFq&>(gf) = polynomial.gf;
	poly = polynomial.poly;

	return *this;
}

// ================================================================================================
GFq_Polynomial& GFq_Polynomial::operator=(const GFq_Element& gfe)
{
	poly.clear();
	const_cast<GFq&>(gf) = gfe.field();
	poly.push_back(gfe);
	return *this;
}

// ================================================================================================
GFq_Polynomial& GFq_Polynomial::operator+=(const GFq_Polynomial& polynomial)
{
	if (gf == polynomial.gf)
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
	}

	return *this;
}

// ================================================================================================
GFq_Polynomial& GFq_Polynomial::operator+=(const GFq_Element& gfe)
{
	poly[0] += gfe;
	return *this;
}

// ================================================================================================
GFq_Polynomial& GFq_Polynomial::operator-=(const GFq_Polynomial& gfe)
{
	return (*this += gfe);
}

// ================================================================================================
GFq_Polynomial& GFq_Polynomial::operator-=(const GFq_Element& gfe)
{
	poly[0] -= gfe;
	return *this;
}

// ================================================================================================
GFq_Polynomial& GFq_Polynomial::operator*=(const GFq_Polynomial& polynomial)
{
	if (gf == polynomial.gf)
	{
		GFq_Polynomial product(gf, deg() + polynomial.deg() + 1);

		for (unsigned int i = 0; i < poly.size(); i++)
		{
			for (unsigned int j = 0; j < polynomial.poly.size(); j++)
			{
				product.poly[i + j] += poly[i] * polynomial.poly[j];
			}
		}

		simplify(product);
		poly = product.poly;
	}

	return *this;
}

// ================================================================================================
GFq_Polynomial& GFq_Polynomial::operator*=(const GFq_Element& gfe)
{
	if (gf == gfe.field())
	{
		for (unsigned int i = 0; i < poly.size(); i++)
		{
			poly[i] *= gfe;
		}
	}

	return *this;
}

// ================================================================================================
GFq_Polynomial& GFq_Polynomial::operator/=(const GFq_Polynomial& divisor)
{
	const std::pair<GFq_Polynomial, GFq_Polynomial>& divres = div(*this, divisor);
	poly = divres.first.get_poly();
	return *this;
}

// ================================================================================================
GFq_Polynomial& GFq_Polynomial::operator/=(const GFq_Element& gfe)
{
	if (gf == gfe.field())
	{
		for (unsigned int i = 0; i < poly.size(); i++)
		{
			poly[i] /= gfe;
		}
	}

	return *this;
}

// ================================================================================================
GFq_Polynomial& GFq_Polynomial::operator%=(const GFq_Polynomial& divisor)
{
	const std::pair<GFq_Polynomial, GFq_Polynomial>& divres = div(*this, divisor);
	poly = divres.second.get_poly();
	return *this;
}

// ================================================================================================
GFq_Polynomial& GFq_Polynomial::operator%=(const unsigned int& power)
{
	if (poly.size() >= power)
	{
		for (size_t i = poly.size(); i > power; i--)
		{
			poly.pop_back();
		}
		//poly.resize(power);
	}

	simplify(*this);
	return *this;
}

// ================================================================================================
GFq_Polynomial& GFq_Polynomial::operator^=(const int& n)
{
    if (n == 0) // P^0 = 1
    {
        poly.clear();
        poly.push_back(GFq_Element(gf, 1)); 
    }
    else if (n > 1)
    {
        GFq_Polynomial result = *this;
        
        for (unsigned int i=0; i<n-1; i++)
        {
            result *= *this;
        }

        *this = result;
    }

	return *this;
}

// ================================================================================================
GFq_Polynomial& GFq_Polynomial::operator<<=(const unsigned int& n)
{
	if (poly.size() > 0)
	{
		std::size_t initial_size = poly.size();
		poly.resize(poly.size() + n, GFq_Element(gf, 1));
		std::copy(poly.rend() - initial_size, poly.rend(), poly.rbegin());
		std::fill(poly.begin(), poly.begin() + n, GFq_Element(gf, 0));
	}

	return *this;
}

// ================================================================================================
GFq_Polynomial& GFq_Polynomial::operator>>=(const unsigned int& n)
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
const GFq_Element& GFq_Polynomial::operator[](const unsigned int& term) const
{
	assert(term < poly.size());
	return poly[term];
}

// ================================================================================================
GFq_Element& GFq_Polynomial::operator[](const unsigned int& term)
{
	assert(term < poly.size());
	return poly[term];
}

// ================================================================================================
GFq_Element GFq_Polynomial::operator()(const GFq_Element& value)
{
	GFq_Element result(gf, 0);

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
const GFq_Element GFq_Polynomial::operator()(const GFq_Element& value) const
{
	GFq_Element result(gf, 0);

	if (poly.size() > 0)
	{
		result = poly[poly.size() - 1];
		for (std::size_t i = poly.size() - 2; static_cast<int>(i) >= 0; i--)
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
GFq_Element GFq_Polynomial::operator()(GFq_Symbol value)
{
	return (*this)(GFq_Element(gf, value));
}

// ================================================================================================
const GFq_Element GFq_Polynomial::operator()(GFq_Symbol value) const
{
	return (*this)(GFq_Element(gf, value));
}

// ================================================================================================
bool GFq_Polynomial::operator==(const GFq_Polynomial& polynomial) const
{
	if (gf == polynomial.gf)
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
	else
		return false;
}

// ================================================================================================
bool GFq_Polynomial::operator!=(const GFq_Polynomial& polynomial) const
{
	return !(*this == polynomial);
}

// ================================================================================================
GFq_Polynomial GFq_Polynomial::derivative()
{
	if ((*this).poly.size() > 1)
	{
		GFq_Polynomial deriv(gf, deg());

		for (unsigned int i = 0; i < poly.size() - 1; i++)
		{
			if (((i + 1) & 1) == 1) // odd power => (bX^n)' = nbX^(n-1) = bX^(n-1)
			{
				deriv.poly[i] = poly[i + 1];
			}
			else
			{
				deriv.poly[i] = 0; // even power => (bX^n)' = nbX^(n-1) = 0X^(n-1) = 0 (eliminated)
			}
		}

		simplify(deriv);
		return deriv;
	}
	else
	{
		return GFq_Polynomial(gf, 0);
	}
}

// ================================================================================================
GFq_Element GFq_Polynomial::make_monic()
{
    GFq_Element lead_coeff = poly.back();
    std::vector<GFq_Element>::iterator it = poly.begin();
    
    for (; it != poly.end(); ++it)
    {
        *it /= lead_coeff;
    }
    
    return lead_coeff;
}

// ================================================================================================
void GFq_Polynomial::rootChien(std::vector<GFq_Element>& roots)
{
	std::vector<GFq_Element> wpoly(poly);
	const GFq_Element zero(gf,0);

	if (poly[0].is_zero())
	{
		roots.push_back(zero);
	}

	for (unsigned int i=0; i < gf.size(); i++)
	{
		GFq_Element sum = std::accumulate(wpoly.begin(), wpoly.end(), zero);

		if (sum.is_zero())
		{
			roots.push_back(GFq_Element(gf,gf.alpha(i)));
		}

		for (unsigned int j=0; j<wpoly.size(); j++)
		{
			wpoly[j] *= gf.alpha(j);
		}
	}
}

// ================================================================================================
void simplify(GFq_Polynomial& polynomial)
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
GFq_Polynomial operator+(const GFq_Polynomial& a, const GFq_Polynomial& b)
{
	GFq_Polynomial result = a;
	result += b;
	return result;
}

// ================================================================================================
GFq_Polynomial operator +(const GFq_Polynomial& a, const GFq_Element& b)
{
	GFq_Polynomial result = a;
	result += b;
	return result;
}

// ================================================================================================
GFq_Polynomial operator +(const GFq_Element& a, const GFq_Polynomial& b)
{
	GFq_Polynomial result = b;
	result += a;
	return result;
}

// ================================================================================================
GFq_Polynomial operator +(const GFq_Polynomial& a, const GFq_Symbol& b)
{
	return a + GFq_Element(a.field(), b);
}

// ================================================================================================
GFq_Polynomial operator +(const GFq_Symbol& a, const GFq_Polynomial& b)
{
	return b + GFq_Element(b.field(), a);
}

// ================================================================================================
GFq_Polynomial operator -(const GFq_Polynomial& a, const GFq_Polynomial& b)
{
	GFq_Polynomial result = a;
	result -= b;
	return result;
}

// ================================================================================================
GFq_Polynomial operator -(const GFq_Polynomial& a, const GFq_Element& b)
{
	GFq_Polynomial result = a;
	result -= b;
	return result;
}

// ================================================================================================
GFq_Polynomial operator -(const GFq_Element& a, const GFq_Polynomial& b)
{
	GFq_Polynomial result = b;
	result -= a;
	return result;
}

// ================================================================================================
GFq_Polynomial operator -(const GFq_Polynomial& a, const GFq_Symbol& b)
{
	return a - GFq_Element(a.field(), b);
}

// ================================================================================================
GFq_Polynomial operator -(const GFq_Symbol& a, const GFq_Polynomial& b)
{
	return b - GFq_Element(b.field(), a);
}

// ================================================================================================
GFq_Polynomial operator *(const GFq_Polynomial& a, const GFq_Polynomial& b)
{
	GFq_Polynomial result = a;
	result *= b;
	return result;
}

// ================================================================================================
GFq_Polynomial operator *(const GFq_Element& a, const GFq_Polynomial& b)
{
	GFq_Polynomial result = b;
	result *= a;
	return result;
}

// ================================================================================================
GFq_Polynomial operator *(const GFq_Polynomial& a, const GFq_Element& b)
{
	GFq_Polynomial result = a;
	result *= b;
	return result;
}

// ================================================================================================
GFq_Polynomial operator /(const GFq_Polynomial& a, const GFq_Polynomial& b)
{
	GFq_Polynomial result = a;
	result /= b;
	return result;
}

// ================================================================================================
GFq_Polynomial operator /(const GFq_Polynomial& a, const GFq_Element& b)
{
	GFq_Polynomial result = a;
	result /= b;
	return result;
}

// ================================================================================================
GFq_Polynomial operator %(const GFq_Polynomial& a, const GFq_Polynomial& b)
{
	GFq_Polynomial result = a;
	result %= b;
	return result;
}

// ================================================================================================
GFq_Polynomial operator %(const GFq_Polynomial& a, const unsigned int& power)
{
	GFq_Polynomial result = a;
	result %= power;
	return result;
}

// ================================================================================================
GFq_Polynomial operator ^(const GFq_Polynomial& a, const int& n)
{
	GFq_Polynomial result = a;
	result ^= n;
	return result;
}

// ================================================================================================
GFq_Polynomial operator<<(const GFq_Polynomial& a, const unsigned int& n)
{
	GFq_Polynomial result = a;
	result <<= n;
	return result;
}

// ================================================================================================
GFq_Polynomial operator>>(const GFq_Polynomial& a, const unsigned int& n)
{
	GFq_Polynomial result = a;
	result >>= n;
	return result;
}

// ================================================================================================
GFq_Polynomial gcd(const GFq_Polynomial& a, const GFq_Polynomial& b)
{
	if ((a.field()) == (b.field()))
	{
		if (a.is_zero() && b.is_zero())
		{
			throw GF_Exception("GCD with both zero operand polynomials");
		}

		if (a.is_zero())
		{
			return b;
		}

		if (b.is_zero())
		{
			return a;
		}

		GFq_Polynomial r = ((a.deg() < b.deg()) ? a : b);
		GFq_Polynomial x = ((a.deg() < b.deg()) ? b : a);
		GFq_Polynomial t = ((a.deg() < b.deg()) ? b : a);

		while (!r.is_zero())
		{
			t = r;
			r = x % t;
			x = t;
		}

		return x;
	}
	else
	{
		throw GF_Exception("GCD with unmatching Galois Fields for operand polynomials");
	}
}

// ================================================================================================
std::pair<GFq_Polynomial, GFq_Polynomial> div(const GFq_Polynomial& dividend,
		const GFq_Polynomial& divisor)
{
	if ((dividend.field() != divisor.field()) || (divisor.deg() < 0))
	{
		throw GF_Exception("GFq Polynomial Division invalid operands");
	}
	else
	{
		if (divisor.deg() == 0)
		{
			GFq_Polynomial quotient(dividend);
			quotient /= divisor[0];
			return std::make_pair(quotient, GFq_Polynomial(quotient.field(), 0));
		}
		else if (dividend.deg() < divisor.deg())
		{
			return std::make_pair(dividend, GFq_Polynomial(dividend.field(), 0));
		}
		else
		{
			GFq_Polynomial remainder(dividend);
			GFq_Polynomial quotient(dividend.field(), dividend.deg() - divisor.deg() + 1);

			while (remainder.is_valid() && (remainder.deg() >= divisor.deg()))
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
std::vector<GFq_Element> rootex_nz(const GFq_Polynomial& a)
{
	std::vector<GFq_Element> roots;
	const GFq& gf = a.field();

	for (unsigned int i=0; i<gf.size(); i++)
	{
		if (a(gf.alpha(i)) == 0)
		{
			roots.push_back(GFq_Element(gf,gf.alpha(i)));
		}
	}

	return roots;
}

// ================================================================================================
std::vector<GFq_Element> rootex(const GFq_Polynomial& a)
{
	std::vector<GFq_Element> roots;
	const GFq& gf = a.field();

	for (unsigned int i=0; i<gf.size()+1; i++)
	{
		if (a(i) == 0)
		{
			roots.push_back(GFq_Element(gf,i));
		}
	}

	return roots;
}

// ================================================================================================
GFq_Polynomial get_monic(const GFq_Polynomial& a, GFq_Element& lead_poly)
{
    GFq_Polynomial result = a;
    lead_poly = result.make_monic();
    return result;
}

// ================================================================================================
std::vector<GFq_Polynomial> square_free_decomposition(const GFq_Polynomial& ff)
{
	const GFq& gf = ff.field();
	GFq_Polynomial f = ff;
	std::vector<GFq_Polynomial> hv;

	unsigned int i = 1;
	GFq_Polynomial u = gcd(f,f.derivative());
	GFq_Polynomial v = f/u;
	GFq_Polynomial w = f.derivative()/u;

	while (!v.is_one())
	{
		GFq_Polynomial h = gcd(v, w-v.derivative());
		v = v/h;
		w = (w-v.derivative())/h;
		hv.push_back(h);
		i++;
	}

	return hv;
}

// ================================================================================================
std::ostream& operator <<(std::ostream& os, const GFq_Polynomial& polynomial)
{
	if (polynomial.deg() >= 0)
	{
		bool is_null = true;
		bool first_coeff = true;

		for (unsigned int i = 0; i < polynomial.poly.size(); i++)
		{
			GFq_Symbol coeff = polynomial.poly[i].poly();

			if (coeff != 0)
			{
				is_null = false;
				os << ((!first_coeff) ? "+ " : "");

				if (polynomial.alpha_format)
				{
					GFq_Symbol log_alpha = polynomial.poly[i].index();

					if (log_alpha == 0)
					{
						if (i == 0)
						{
							os << "1 ";
						}
						else
						{
							os << "";
						}
					}
					else if (log_alpha == 1)
					{
						os << "a" << ((!first_coeff) ? "*" : " ");
					}
					else
					{
						os << "a^" << log_alpha << ((!first_coeff) ? "*" : " ");
					}
				}
				else
				{
					if (coeff == 1)
					{
						if (i != 0)
						{
							os << "";
						}
						else
						{
							os << coeff << " ";
						}
					}
					else
					{
						os << coeff << ((!first_coeff) ? "*" : " ");
					}
				}

				if (i != 0)
				{
					if (i == 1)
					{
						os << "X ";
					}
					else
					{
						os << "X^" << i << " ";
					}
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
