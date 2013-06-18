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
#include "GF_Utils.h"
#include <set>

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
GFq_BivariatePolynomial::GFq_BivariatePolynomial(const std::pair<unsigned int, unsigned int>& _weights) :
        weights(_weights),
        monomials(GFq_WeightedRevLex_BivariateMonomial(_weights))
{}

// ================================================================================================
GFq_BivariatePolynomial::GFq_BivariatePolynomial(const GFq_BivariatePolynomial& polynomial) :
		weights(polynomial.get_weights()),
		monomials(GFq_WeightedRevLex_BivariateMonomial(polynomial.get_weights()))
{
    monomials = polynomial.monomials;
}

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
    monomials = polynomial.monomials;
}

// ================================================================================================
void GFq_BivariatePolynomial::init(const std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& _monomials)
{
	monomials = _monomials;
}

// ================================================================================================
void GFq_BivariatePolynomial::init_x_pow(const GFq& gf, unsigned int x_pow)
{
    monomials.insert(std::make_pair(std::make_pair(x_pow,0), GFq_Element(gf, 1)));
}

// ================================================================================================
void GFq_BivariatePolynomial::init_y_pow(const GFq& gf, unsigned int y_pow)
{
    monomials.insert(std::make_pair(std::make_pair(0, y_pow), GFq_Element(gf, 1)));
}

// ================================================================================================
void GFq_BivariatePolynomial::init_x_pow_series(const GFq& gf, unsigned int max_pow)
{
    for (unsigned int i = 0; i <= max_pow; i++)
    {
        monomials.insert(std::make_pair(std::make_pair(i,0), GFq_Element(gf, 1)));
    }
}

// ================================================================================================
void GFq_BivariatePolynomial::init_y_pow_series(const GFq& gf, unsigned int max_pow)
{
    for (unsigned int i = 0; i <= max_pow; i++)
    {
        monomials.insert(std::make_pair(std::make_pair(0,i), GFq_Element(gf, 1)));
    }
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
unsigned int GFq_BivariatePolynomial::wdeg() const
{
    std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::const_iterator mono_it = monomials.begin();
    unsigned int max_wd = 0;
    
    for (; mono_it != monomials.end(); ++mono_it)
    {
        GFq_BivariateMonomial mono = static_cast<GFq_BivariateMonomial>(*mono_it);
        
        if (mono.wdeg(weights) > max_wd)
        {
            max_wd = mono.wdeg(weights);
        }
    }
    
    return max_wd;
}

// ================================================================================================
void GFq_BivariatePolynomial::sum(std::vector<GFq_BivariateMonomial>& _sum_monomials, const GFq_BivariatePolynomial& a, const GFq_BivariatePolynomial& b)
{
	if (a.get_weights() != b.get_weights())
	{
		throw GF_Exception("Cannot add bivariate polynomials with different degree weights");
	}
	else
	{
		std::vector<GFq_BivariateMonomial> sum_monomials;
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

        // simplify (removes monomials with zero coefficients)

		std::vector<GFq_BivariateMonomial>::iterator mono_it = sum_monomials.begin();

        for (; mono_it != sum_monomials.end(); ++mono_it)
        {
            if (!(mono_it->second.is_zero()))
            {
            	_sum_monomials.push_back(*mono_it);
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
		product(prod_monomials, a.get_monomials(), b.get_monomials());
		simplify(prod_monomials); // simplify (removes monomials with zero coefficients)
	}
}

// ================================================================================================
void GFq_BivariatePolynomial::product(std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& prod_monomials,
		const std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& a_monomials,
		const std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& b_monomials)
{
	std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::const_iterator a_it = a_monomials.begin();

	for (; a_it != a_monomials.end(); ++a_it)
	{
		std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::const_iterator b_it = b_monomials.begin();

		for (; b_it != b_monomials.end(); ++b_it)
		{
			GFq_BivariateMonomial mono_product = static_cast<GFq_BivariateMonomial>(*a_it) * static_cast<GFq_BivariateMonomial>(*b_it);
			add_monomial(prod_monomials, mono_product);
		}
	}
}

// ================================================================================================
void GFq_BivariatePolynomial::product(std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& prod_monomials,
		const GFq_Element& v,
		const std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& a_monomials,
		const std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& b_monomials)
{
	std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::const_iterator a_it = a_monomials.begin();

	for (; a_it != a_monomials.end(); ++a_it)
	{
		std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::const_iterator b_it = b_monomials.begin();

		for (; b_it != b_monomials.end(); ++b_it)
		{
			GFq_BivariateMonomial mono_product = v * static_cast<GFq_BivariateMonomial>(*a_it) * static_cast<GFq_BivariateMonomial>(*b_it);
			add_monomial(prod_monomials, mono_product);
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
void GFq_BivariatePolynomial::pow(std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& pow_monomials,
		const GFq_BivariatePolynomial& a,
		unsigned int n)
{
	if (!a.is_valid())
	{
		throw GF_Exception("Invalid polynomial");
	}
	else
	{
		const std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& a_monomials = a.get_monomials();
		const GFq& gf = a_monomials.begin()->second.field();

		if (n == 0)
		{
			pow_monomials.insert(std::make_pair(std::make_pair(0,0), GFq_Element(gf,1))); // P^n = 1
		}
		else if (n == 1) // P^1 = P
		{
			pow_monomials = a_monomials;
		}
		else if (a_monomials.size() == 1) // Optimization if polynomial is reduced to a single monomial: (a*X^i*Y^j)^n = a^n*X^(i*n)*Y^(j*n)
        {
            add_monomial(pow_monomials, (a_monomials.begin()->second)^n, (a_monomials.begin()->first.first)*n, (a_monomials.begin()->first.second)*n);
        }
        else
		{
			std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial> at_monomials(a_monomials); //!< temporary a^i monomials
			std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial> pt_monomials(pow_monomials); //!< temporary power monomials

			for (unsigned int i=0; i<n-1; i++)
			{
				pt_monomials.clear();

				std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::const_iterator at_it = at_monomials.begin();
				for (; at_it != at_monomials.end(); ++at_it)
				{
					std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::const_iterator a_it = a_monomials.begin();
					for (; a_it != a_monomials.end(); ++a_it)
					{
						GFq_BivariateMonomial mono_product = static_cast<GFq_BivariateMonomial>(*at_it) * static_cast<GFq_BivariateMonomial>(*a_it);
						std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::iterator mono_it = pt_monomials.find(mono_product.get_exponents());

						if (mono_it == pt_monomials.end()) // exponents pair does not exist yet
						{
							//std::cout << "new: " << *at_it << " * " << *a_it << " = " << mono_product << std::endl;
							pt_monomials.insert(mono_product); // insert the new monomial
						}
						else // exponents pair already exist
						{
							//std::cout << "old: " << *at_it << " * " << *a_it << " = " << mono_product << std::endl;
							mono_it->second += mono_product.second; // add coefficients
							//std::cout << " -> " << *mono_it << std::endl;
						}
					}
				}

				at_monomials = pt_monomials;
			}

			// copy last temporary a^i monomials to result and simplify

	        std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::iterator mono_it = at_monomials.begin();

	        for (; mono_it != at_monomials.end(); ++mono_it)
	        {
	            if (!mono_it->second.is_zero())
	            {
	            	add_monomial(pow_monomials, *mono_it);
	            }
	        }
		}
	}
}

// ================================================================================================
void GFq_BivariatePolynomial::add_monomial(std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& monomials,
		const GFq_BivariateMonomialKeyValueRepresentation& a)
{
	std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::iterator mono_it = monomials.find(a.first);

	if (mono_it == monomials.end()) // exponents pair does not exist yet
	{
		monomials.insert(a); // insert the new monomial
	}
	else // exponents pair already exist
	{
		mono_it->second += a.second; // add coefficients
	}
}

// ================================================================================================
void GFq_BivariatePolynomial::add_monomial(std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& monomials,
		const GFq_Element& coeff,
		unsigned int x_pow,
		unsigned int y_pow)
{
	GFq_BivariateMonomialExponents exponents(x_pow, y_pow);
	std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::iterator mono_it = monomials.find(exponents);

	if (mono_it == monomials.end()) // exponents pair does not exist yet
	{
		GFq_BivariateMonomialKeyValueRepresentation a(exponents, coeff);
		monomials.insert(a); // insert the new monomial
	}
	else // exponents pair already exist
	{
		mono_it->second += coeff; // add coefficients
	}
}

// ================================================================================================
void GFq_BivariatePolynomial::simplify(std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& monomials)
{
    std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::iterator mono_it = monomials.begin();

    while (mono_it != monomials.end())
    {
        if (mono_it->second.is_zero())
        {
        	monomials.erase(mono_it++);
        }
        else
        {
        	++mono_it;
        }
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
GFq_BivariatePolynomial& GFq_BivariatePolynomial::operator*=(const GFq_BivariateMonomial& monomial)
{
	GFq_WeightedRevLex_BivariateMonomial mono_exp_compare(weights);
	std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial> product_monomials(mono_exp_compare);
	std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::const_iterator mono_it = monomials.begin();
	
	for (; mono_it != monomials.end(); ++mono_it)
	{
		unsigned int eX = mono_it->first.first + monomial.eX();
		unsigned int eY = mono_it->first.second + monomial.eY();
		product_monomials.insert(std::make_pair(std::make_pair(eX, eY), mono_it->second * monomial.coeff()));
	}
	
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
GFq_BivariatePolynomial& GFq_BivariatePolynomial::operator^=(unsigned int n)
{
	GFq_WeightedRevLex_BivariateMonomial mono_exp_compare(weights);
	std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial> pow_monomials(mono_exp_compare);

	pow(pow_monomials, *this, n);
	monomials = pow_monomials;

	return *this;
}


// ================================================================================================
bool GFq_BivariatePolynomial::operator==(const GFq_BivariatePolynomial& polynomial) const
{
	return (weights == polynomial.weights) && (monomials == polynomial.monomials);
}

// ================================================================================================
bool GFq_BivariatePolynomial::operator!=(const GFq_BivariatePolynomial& polynomial) const
{
	return (weights != polynomial.weights) || (monomials != polynomial.monomials);
}

// ================================================================================================
GFq_Element GFq_BivariatePolynomial::operator()(const GFq_Element& x_value, const GFq_Element& y_value) const
{
	if (x_value.field() != y_value.field())
	{
		throw GF_Exception("point coordinates must be of the same Galois Field to evaluate bivariate polynomial at this point");
	}
	else
	{
		std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::const_iterator mono_it = monomials.begin();
		GFq_Element result(mono_it->second.field(), 0);

		for (; mono_it != monomials.end(); ++mono_it)
		{
			result += (x_value^mono_it->first.first) * (y_value^mono_it->first.second) * (mono_it->second);
		}

		return result;
	}
}

// ================================================================================================
GFq_Polynomial GFq_BivariatePolynomial::get_X_0() const
{
	return get_v_0(true);
}

// ================================================================================================
GFq_Polynomial GFq_BivariatePolynomial::get_0_Y() const
{
	return get_v_0(false);
}

// ================================================================================================
GFq_Polynomial GFq_BivariatePolynomial::get_v_0(bool x_terms) const
{
	if (monomials.size() == 0)
	{
		throw GF_Exception("Bivariate polynomial is invalid");
	}
	else
	{
		std::map<unsigned int, GFq_Element> poly_map;
		std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::const_iterator mono_it = monomials.begin();
		const GFq& gf = mono_it->second.field();
		GFq_Element zero(gf,0);

		for (; mono_it != monomials.end(); ++mono_it)
		{
			if (x_terms)
			{
				if (mono_it->first.second == 0) // no Y term
				{
					poly_map.insert(std::make_pair(mono_it->first.first, mono_it->second));
				}
			}
			else
			{
				if (mono_it->first.first == 0) // no X term
				{
					poly_map.insert(std::make_pair(mono_it->first.second, mono_it->second));
				}
			}
		}

		if (poly_map.size() == 0)
		{
			return GFq_Polynomial(zero);
		}
		else
		{
			unsigned int max_pow = poly_map.rbegin()->first;
			std::vector<GFq_Element> poly(max_pow+1, zero);
			std::map<unsigned int, GFq_Element>::const_iterator x_map_it = poly_map.begin();

			for (; x_map_it != poly_map.end(); ++x_map_it)
			{
				poly[x_map_it->first] = x_map_it->second;
			}

			return GFq_Polynomial(gf, poly);
		}
	}
}

// ================================================================================================
GFq_BivariatePolynomial GFq_BivariatePolynomial::operator()(const GFq_BivariatePolynomial& P, const GFq_BivariatePolynomial& Q) const
{
	if (monomials.size() == 0)
	{
		throw GF_Exception("Bivariate polynomial is invalid");
	}
	else if (!P.is_valid())
	{
		throw GF_Exception("First operand polynomial is invalid");
	}
	else if (!Q.is_valid())
	{
		throw GF_Exception("Second operand polynomial is invalid");
	}
	else if (P.get_weights() != Q.get_weights())
	{
		throw GF_Exception("Cannot evaluate with polynomials of different weights");
	}
	else
	{
		std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::const_iterator mono_it = monomials.begin();
		const GFq& gf = mono_it->second.field();
		GFq_WeightedRevLex_BivariateMonomial mono_exp_compare(weights); // Use reverse lexical order
		std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial> result_monomials(mono_exp_compare);

		for (; mono_it != monomials.end(); ++mono_it)
		{
			std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial> partX_monomials(mono_exp_compare);
			pow(partX_monomials, P, mono_it->first.first); // X^i -> (P(X,Y))^i
			std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial> partY_monomials(mono_exp_compare);
			pow(partY_monomials, Q, mono_it->first.second); // Y^j -> (Q(X,Y))^j
			product(result_monomials, mono_it->second, partX_monomials, partY_monomials); // a_i,j*(P(X,Y))^i*(Q(X,Y))^j
		}

		simplify(result_monomials); // simplify the resulting monomials (removes those with zero coefficient)
		GFq_BivariatePolynomial result(weights);
		result.init(result_monomials);
		return result;
	}
}

// ================================================================================================
GFq_BivariatePolynomial& GFq_BivariatePolynomial::make_star()
{
	if (monomials.size() == 0)
	{
		throw GF_Exception("Bivariate polynomial is invalid");
	}
	else
	{
		std::set<unsigned int> x_exp_set;
		std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::const_iterator mono_it = monomials.begin();
		const GFq& gf = mono_it->second.field();
		
		for (; mono_it != monomials.end(); ++mono_it)
		{
			x_exp_set.insert(mono_it->first.first);
		}
		
		unsigned int h = *(x_exp_set.begin());
		
		if (h > 0)
		{
			GFq_WeightedRevLex_BivariateMonomial mono_exp_compare(weights);
			std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial> div_monomials(mono_exp_compare);
			mono_it = monomials.begin();
			
			for (; mono_it != monomials.end(); ++mono_it)
			{
				div_monomials.insert(std::make_pair(std::make_pair(mono_it->first.first-h, mono_it->first.second), mono_it->second));
			}
			
			monomials = div_monomials;
		}
	}
	
	return *this;
}

// ================================================================================================
GFq_BivariatePolynomial& GFq_BivariatePolynomial::make_dHasse(unsigned int mu, unsigned int nu)
{
	if (monomials.size() == 0)
	{
		throw GF_Exception("Bivariate polynomial is invalid");
	}
	else if ((mu != 0) or (nu != 0)) // ^[0,0] is the trivial case where the polynomial is unmodified
	{
		GFq_WeightedRevLex_BivariateMonomial mono_exp_compare(weights);
		std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial> hasse_monomials(mono_exp_compare);
		std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>::const_iterator mono_it = monomials.begin();
		const GFq& gf = mono_it->second.field();

		for (; mono_it != monomials.end(); ++mono_it)
		{
			unsigned int eX = mono_it->first.first;
			unsigned int eY = mono_it->first.second;

			if (!(eX < mu) && !(eY < nu))
			{
				// coefficient is not null if both binomial coefficients are even
				if (!(binomial_coeff_parity(eX,mu)) && !(binomial_coeff_parity(eY,nu)))
				{
					hasse_monomials.insert(std::make_pair(std::make_pair(eX-mu, eY-nu), mono_it->second));
				}
			}
		}

		if (hasse_monomials.size() == 0)
		{
			monomials.clear();
			monomials.insert(std::make_pair(std::make_pair(0,0), GFq_Element(gf, 0)));
		}
		else
		{
			monomials = hasse_monomials;
		}
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
	std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial>& monomials = polynomial.get_monomials_for_update();
	GFq_BivariatePolynomial::simplify(monomials);
}

// ================================================================================================
GFq_BivariatePolynomial operator+(const GFq_BivariatePolynomial& a, const GFq_BivariatePolynomial& b)
{
	GFq_BivariatePolynomial result(a.get_weights()); // copy without creating monomials
	std::vector<GFq_BivariateMonomial> sum_monomials;

	GFq_BivariatePolynomial::sum(sum_monomials, a, b);

	result.init(sum_monomials);

	return result;
}

// ================================================================================================
GFq_BivariatePolynomial operator +(const GFq_BivariatePolynomial& a, const GFq_Element& b)
{
	GFq_BivariatePolynomial result(a);
	result += b;
	return result;
}

// ================================================================================================
GFq_BivariatePolynomial operator +(const GFq_Element& a, const GFq_BivariatePolynomial& b)
{
	GFq_BivariatePolynomial result(b);
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
	GFq_BivariatePolynomial result(a.get_weights()); // copy without creating monomials
	GFq_BivariatePolynomial::product(product_monomials, a, b);
	result.init(product_monomials);

	return result;
}

// ================================================================================================
GFq_BivariatePolynomial operator*(const GFq_BivariatePolynomial& a, const GFq_BivariateMonomial& b)
{
	GFq_BivariatePolynomial result(a);
	result *= b;
	return result;
}

// ================================================================================================
GFq_BivariatePolynomial operator*(const GFq_BivariatePolynomial& a, const GFq_Element& b)
{
	GFq_BivariatePolynomial result(a);
	result *= b;
	return result;
}

// ================================================================================================
GFq_BivariatePolynomial operator*(const GFq_Element& a, const GFq_BivariatePolynomial& b)
{
	GFq_BivariatePolynomial result(b);
	result *= a;
	return result;
}

// ================================================================================================
GFq_BivariatePolynomial operator /(const GFq_BivariatePolynomial& a, const GFq_BivariateMonomial& b)
{
	GFq_WeightedRevLex_BivariateMonomial mono_exp_compare(a.get_weights());
	std::map<GFq_BivariateMonomialExponents, GFq_Element, GFq_WeightedRevLex_BivariateMonomial> div_monomials(mono_exp_compare);
	GFq_BivariatePolynomial result(a.get_weights()); // copy without creating monomials
	GFq_BivariatePolynomial::division(div_monomials, a, b);
	result.init(div_monomials);

	return result;
}

// ================================================================================================
GFq_BivariatePolynomial operator /(const GFq_BivariatePolynomial& a, const GFq_Element& b)
{
	GFq_BivariatePolynomial result(a);
	result /= b;
	return result;
}

// ================================================================================================
GFq_BivariatePolynomial operator ^(const GFq_BivariatePolynomial& a, unsigned int n)
{
	GFq_BivariatePolynomial result(a);
	result ^= n;
	return result;
}

// ================================================================================================
GFq_BivariatePolynomial star(const GFq_BivariatePolynomial& a)
{
	GFq_BivariatePolynomial result(a);
	result.make_star();
	return result;
}

// ================================================================================================
GFq_BivariatePolynomial dHasse(unsigned int mu, unsigned int nu, const GFq_BivariatePolynomial& a)
{
	GFq_BivariatePolynomial result(a);
	result.make_dHasse(mu, nu);
	return result;
}

} // namespace gf
} // namespace rssoft
