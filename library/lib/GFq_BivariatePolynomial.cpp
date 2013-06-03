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

 Bivariate polynomials with coefficients in GF(2^m) class

 */

#include "GFq_BivariatePolynomial.h"

namespace rssoft
{
namespace gf
{

// ================================================================================================
GFq_BivariatePolynomial::GFq_BivariatePolynomial(unsigned int w_x, unsigned int w_y) :
		weights(w_x,w_y),
		monomials(GFq_WeightedRevLex_BivariateMonomial(w_x,w_y))
{}

// ================================================================================================
GFq_BivariatePolynomial::GFq_BivariatePolynomial(const GFq_BivariatePolynomial& polynomial) :
		weights(polynomial.get_weights()),
		monomials(GFq_WeightedRevLex_BivariateMonomial(polynomial.get_weights()))
{}

// ================================================================================================
GFq_BivariatePolynomial::~GFq_BivariatePolynomial()
{}

// ================================================================================================
void GFq_BivariatePolynomial::init(std::vector<GFq_BivariateMonomial>& _monomials)
{
	monomials.clear();
	monomials.insert(_monomials.begin(), _monomials.end());
}

// ================================================================================================
void GFq_BivariatePolynomial::init(const GFq_BivariatePolynomial& polynomial)
{
	monomials.clear();
	monomials.insert(polynomial.monomials.begin(), polynomial.monomials.end());
}

// ================================================================================================
void GFq_BivariatePolynomial::init(const std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& _monomials)
{
	monomials = _monomials;
}

// ================================================================================================
bool GFq_BivariatePolynomial::is_valid() const
{
	return monomials.size() != 0;
}

// ================================================================================================
bool GFq_BivariatePolynomial::is_const(GFq_Element& const_value) const
{
	if (monomials.size() == 1)
	{
		const GFq_BivariateMonomialKeyValueRepresentation& m0 = *(monomials.begin());

		if (m0.second == const_value)
		{
			return (m0.first.are_zero());
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
bool GFq_BivariatePolynomial::is_zero() const
{
	if (monomials.size() == 0)
	{
		return true;
	}
	else
	{
		if (monomials.size() == 1)
		{
			const GFq_BivariateMonomialKeyValueRepresentation& m0 = *(monomials.begin());

			if (m0.second.is_zero())
			{
				return (m0.first.are_zero());
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
}

// ================================================================================================
bool GFq_BivariatePolynomial::is_one() const
{
	if (monomials.size() == 0)
	{
		return true;
	}
	else
	{
		if (monomials.size() == 1)
		{
			const GFq_BivariateMonomialKeyValueRepresentation& m0 = *(monomials.begin());

			if (m0.second.is_one())
			{
				return (m0.first.are_zero());
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
}

// ================================================================================================
GFq_BivariatePolynomial& GFq_BivariatePolynomial::operator =(const GFq_BivariatePolynomial& polynomial)
{
	weights = polynomial.weights;
	monomials = polynomial.monomials;
	return *this;
}


// ================================================================================================
void GFq_BivariatePolynomial::sum(std::vector<GFq_BivariateMonomial>& sum_monomials, const GFq_BivariatePolynomial& a, const GFq_BivariatePolynomial& b)
{
	if (a.get_weights() != b.get_weights())
	{
		throw GF_Exception("Cannot add bivariate polynomials with different degree weights");
	}
	else
	{
		GFq_WeightedRevLex_BivariateMonomial mono_exp_compare(a.get_weights()); // Use reverse lexical order
		const std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& a_monomials = a.get_monomials();
		const std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& b_monomials = b.get_monomials();
		std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::const_iterator a_it = a_monomials.begin();
		std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::const_iterator b_it = b_monomials.begin();

		while ((a_it != a_monomials.end()) || (b_it != b_monomials.end()))
		{
			if (a_it == a_monomials.end())
			{
				//std::cout << "tail copy b " << *b_it << std::endl;
				sum_monomials.push_back(static_cast<GFq_BivariateMonomial>(*b_it));
				++b_it;
			}
			else if (b_it ==  b_monomials.end())
			{
				//std::cout << "tail copy a " << *a_it << std::endl;
				sum_monomials.push_back(static_cast<GFq_BivariateMonomial>(*a_it));
				++a_it;
			}
			else
			{
				if (mono_exp_compare(a_it->first, b_it->first)) // a monomial is less than b
				{
					//std::cout << "copy a " << *a_it << " b is " << *b_it << std::endl;
					sum_monomials.push_back(static_cast<GFq_BivariateMonomial>(*a_it));
					++a_it;
				}
				else if (mono_exp_compare(b_it->first, a_it->first)) // b monomial is less than a
				{
					//std::cout << "copy b " << *b_it << " a is " << *a_it << std::endl;
					sum_monomials.push_back(static_cast<GFq_BivariateMonomial>(*b_it));
					++b_it;
				}
				else // monomial orders are equal (thus their exponents are equal) so coefficient can be summed up
				{
					//std::cout << "copy sum " << *a_it << " + " << *b_it << std::endl;
					sum_monomials.push_back(static_cast<GFq_BivariateMonomial>(*a_it));
					sum_monomials.back().second += b_it->second;
					++a_it;
					++b_it;
				}
			}
		}
	}
}


// ================================================================================================
void GFq_BivariatePolynomial::product(std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& prod_monomials,
			const GFq_BivariatePolynomial& a,
			const GFq_BivariatePolynomial& b)
{
	if (a.get_weights() != b.get_weights())
	{
		throw GF_Exception("Cannot add bivariate polynomials with different degree weights");
	}
	else
	{
		const std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& a_monomials = a.get_monomials();
		const std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& b_monomials = b.get_monomials();
		std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::const_iterator a_it = a_monomials.begin();

		for (; a_it != a_monomials.end(); ++a_it)
		{
			std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::const_iterator b_it = b_monomials.begin();

			for (; b_it != b_monomials.end(); ++b_it)
			{
				GFq_BivariateMonomial mono_product = static_cast<GFq_BivariateMonomial>(*a_it) * static_cast<GFq_BivariateMonomial>(*b_it);
				std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::iterator mono_it = prod_monomials.find(mono_product.get_exponents());

				if (mono_it == prod_monomials.end()) // exponents pair does not exist yet
				{
					prod_monomials.insert(mono_product); // insert the new monomial
				}
				else // exponents pair already exist
				{
					mono_it->second += mono_product.second; // add coefficients
				}

			}
		}
	}
}

// ================================================================================================
void GFq_BivariatePolynomial::division(std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& div_monomials,
		const GFq_BivariatePolynomial& a,
		const GFq_BivariateMonomial& b)
{
	const std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& a_monomials = a.get_monomials();
	std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::const_iterator a_it = a_monomials.begin();

	for (; a_it != a_monomials.end(); ++a_it)
	{
		GFq_BivariateMonomial mono_div = static_cast<GFq_BivariateMonomial>(*a_it) / b;
		div_monomials.insert(mono_div);
	}
}

// ================================================================================================
GFq_BivariatePolynomial& GFq_BivariatePolynomial::operator+=(const GFq_BivariatePolynomial& polynomial)
{
	std::vector<GFq_BivariateMonomial> sum_monomials;

	sum(sum_monomials, *this, polynomial);

	monomials.clear();
	monomials.insert(sum_monomials.begin(), sum_monomials.end());

	return *this;
}

// ================================================================================================
GFq_BivariatePolynomial& GFq_BivariatePolynomial::operator+=(const GFq_Element& gfe)
{
	if (monomials.size() == 0)
	{
		monomials.insert(std::make_pair(std::make_pair(0,0), gfe));
	}
	else
	{
		if (monomials.begin()->first.are_zero()) // constant monomial
		{
			monomials.begin()->second += gfe;
		}
		else
		{
			monomials.insert(std::make_pair(std::make_pair(0,0), gfe));
		}
	}

	return *this;
}

// ================================================================================================
GFq_BivariatePolynomial& GFq_BivariatePolynomial::operator-=(const GFq_BivariatePolynomial& polynomial)
{
	return (*this += polynomial);
}

// ================================================================================================
GFq_BivariatePolynomial& GFq_BivariatePolynomial::operator-=(const GFq_Element& gfe)
{
	return (*this += gfe);
}

// ================================================================================================
GFq_BivariatePolynomial& GFq_BivariatePolynomial::operator*=(const GFq_BivariatePolynomial& polynomial)
{
	GFq_WeightedRevLex_BivariateMonomial mono_exp_compare(weights);
	std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial> product_monomials(mono_exp_compare);

	product(product_monomials, *this, polynomial);
	monomials = product_monomials;

	return *this;
}

// ================================================================================================
GFq_BivariatePolynomial& GFq_BivariatePolynomial::operator*=(const GFq_Element& gfe)
{
	std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::iterator mono_it = monomials.begin();

	for (; mono_it != monomials.end(); ++mono_it)
	{
		mono_it->second *= gfe;
	}

	return *this;
}

// ================================================================================================
GFq_BivariatePolynomial& GFq_BivariatePolynomial::operator/=(const GFq_BivariateMonomial& monomial)
{
	GFq_WeightedRevLex_BivariateMonomial mono_exp_compare(weights);
	std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial> div_monomials(mono_exp_compare);

	division(div_monomials, *this, monomial);
	monomials = div_monomials;

	return *this;
}

// ================================================================================================
GFq_BivariatePolynomial& GFq_BivariatePolynomial::operator/=(const GFq_Element& gfe)
{
	std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::iterator mono_it = monomials.begin();

	for (; mono_it != monomials.end(); ++mono_it)
	{
		mono_it->second /= gfe;
	}

	return *this;
}

// ================================================================================================
std::ostream& operator <<(std::ostream& os, const GFq_BivariatePolynomial& polynomial)
{
	if (polynomial.monomials.size() == 0)
	{
		os << "<invalid>";
	}
	else
	{
		std::map<GFq_BivariateMonomialExponents, GFq_Element>::const_iterator it = polynomial.monomials.begin();

		for (; it != polynomial.monomials.end(); ++it)
		{
			if (it != polynomial.monomials.begin())
			{
				os << " + ";
			}

			os << *it;
		}
	}

	return os;
}

// ================================================================================================
void simplify(GFq_BivariatePolynomial& polynomial)
{
	const std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& monomials = polynomial.get_monomials();
	std::vector<GFq_BivariateMonomial> sum_monomials;

	if (monomials.size() > 0)
	{
		std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::const_iterator it = monomials.begin();

		for (; it != monomials.end(); ++it)
		{
			if (!(it->second.is_zero()))
			{
				sum_monomials.push_back(static_cast<GFq_BivariateMonomial>(*it));
			}
		}
	}

	polynomial.init(sum_monomials);
}

// ================================================================================================
GFq_BivariatePolynomial operator+(const GFq_BivariatePolynomial& a, const GFq_BivariatePolynomial& b)
{
	GFq_BivariatePolynomial result(a);
	std::vector<GFq_BivariateMonomial> sum_monomials;

	GFq_BivariatePolynomial::sum(sum_monomials, a, b);

	result.init(sum_monomials);

	return result;
}

// ================================================================================================
GFq_BivariatePolynomial operator +(const GFq_BivariatePolynomial& a, const GFq_Element& b)
{
	GFq_BivariatePolynomial result(a);
	result.init(a);
	result += b;
	return result;
}

// ================================================================================================
GFq_BivariatePolynomial operator +(const GFq_Element& a, const GFq_BivariatePolynomial& b)
{
	GFq_BivariatePolynomial result(b);
	result.init(b);
	result += a;
	return result;
}

// ================================================================================================
GFq_BivariatePolynomial operator-(const GFq_BivariatePolynomial& a, const GFq_BivariatePolynomial& b)
{
	return a+b;
}

// ================================================================================================
GFq_BivariatePolynomial operator -(const GFq_BivariatePolynomial& a, const GFq_Element& b)
{
	return a+b;
}

// ================================================================================================
GFq_BivariatePolynomial operator -(const GFq_Element& a, const GFq_BivariatePolynomial& b)
{
	return a+b;
}

// ================================================================================================
GFq_BivariatePolynomial operator*(const GFq_BivariatePolynomial& a, const GFq_BivariatePolynomial& b)
{
	GFq_WeightedRevLex_BivariateMonomial mono_exp_compare(a.get_weights());
	std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial> product_monomials(mono_exp_compare);
	GFq_BivariatePolynomial result(a);
	GFq_BivariatePolynomial::product(product_monomials, a, b);
	result.init(product_monomials);

	return result;
}

// ================================================================================================
GFq_BivariatePolynomial operator*(const GFq_BivariatePolynomial& a, const GFq_Element& b)
{
	GFq_BivariatePolynomial result(a);
	result.init(a);
	result *= b;
	return result;
}

// ================================================================================================
GFq_BivariatePolynomial operator*(const GFq_Element& a, const GFq_BivariatePolynomial& b)
{
	GFq_BivariatePolynomial result(b);
	result.init(b);
	result *= a;
	return result;
}

// ================================================================================================
GFq_BivariatePolynomial operator /(const GFq_BivariatePolynomial& a, const GFq_BivariateMonomial& b)
{
	GFq_WeightedRevLex_BivariateMonomial mono_exp_compare(a.get_weights());
	std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial> div_monomials(mono_exp_compare);
	GFq_BivariatePolynomial result(a);
	GFq_BivariatePolynomial::division(div_monomials, a, b);
	result.init(div_monomials);

	return result;
}

// ================================================================================================
GFq_BivariatePolynomial operator /(const GFq_BivariatePolynomial& a, const GFq_Element& b)
{
	GFq_BivariatePolynomial result(a);
	result.init(a);
	result /= b;
	return result;
}

} // namespace gf
} // namespace rssoft
