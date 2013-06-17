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

 Roth-Ruckenstein factorization class for soft decision decoding
 Optimized recursive strategy

 */
#ifndef __RR_FACTORIZATION_H__
#define __RR_FACTORIZATION_H__

#include "GFq_Polynomial.h"
#include "GFq_Element.h"
#include <vector>
#include <set>

namespace rssoft
{

namespace gf
{
class GFq;
class GFq_BivariatePolynomial;
}

/**
 * \brief Node in the Roth-Ruckenstein's algorithm
 */
class RR_Node
{
public:
	/**
	 * Constructor
	 * \param _parent Pointer to the parent node, 0 for root node
	 * \param _Q Node polynomial
     * \param _coeff Coefficient on the arc towards this node
	 * \param _id Node identifier
	 */
	RR_Node(RR_Node *_parent,
			const gf::GFq_BivariatePolynomial& _Q,
            const gf::GFq_Element& _coeff,
			unsigned int _id);

	/**
	 * Get the node's Id
	 */
	unsigned int get_id() const
	{
		return id;
	}

	/**
	 * Get the degree of the node
	 */
	int get_degree() const
	{
		return degree;
	}

	/**
	 * Get node's polynomial
	 */
	const gf::GFq_BivariatePolynomial& getQ() const
	{
		return Q;
	}

	/**
	 * Get coefficient towards the node
	 */
	const gf::GFq_Element& get_coeff() const
	{
		return coeff;
	}

	/**
	 * Add a root in Y
	 */
	void add_ry(const gf::GFq_Element& root_y)
	{
		ry_set.insert(root_y);
	}

	/**
	 * Tells if an element is in the set of Y roots
	 */
	bool is_in_ry_set(const gf::GFq_Element& root_y) const
	{
		std::set<gf::GFq_Element>::const_iterator e_it = ry_set.find(root_y);
		return e_it != ry_set.end();
	}

protected:
	RR_Node *parent; //!< Pointer to the parent node
	const gf::GFq_BivariatePolynomial& Q; //!< Node's polynomial
    const gf::GFq_Element& coeff; // !< Coefficient on the arc towards this node
	unsigned int id; //!< Identifier number of the node
	int degree; //!< The distance of the node from the root counted in the number of arcs less one
	std::set<gf::GFq_Element> ry_set; //!< Set of roots in Y of the node's polynomial

};

/**
 * \brief Roth-Ruckenstein's factorization
 */
class RR_Factorization
{
public:
	/**
	 * Constructor
	 * \param gf Reference to the Galois Field being used
	 * \param k as in RS(n,k)
	 */
	RR_Factorization(const gf::GFq& _gf, unsigned int _k);

	/**
	 * Destructor
	 */
	~RR_Factorization();

	/**
	 * Initialization before run
	 */
	void init();

	/**
	 * Run factorization of given polynomial
	 * \param polynomial Input polynomial
     * \return list of polynomial factors
	 */
	std::vector<gf::GFq_Polynomial>& run(const gf::GFq_BivariatePolynomial& polynomial);

protected:
	/**
	 * Recursive run on a node
	 * \param rr_node The node
	 * \return partial polynomial built at the node
	 */
	gf::GFq_Polynomial node_run(RR_Node& rr_node);

	const gf::GFq& gf; //!< Reference to the Galois Field being used
	unsigned int k;    //!< k as in RS(n,k)
	unsigned int t;    //!< nodes but root node count
	std::vector<gf::GFq_Polynomial> F; //!< Result list of f(X) polynomials
};

} // namespace rssoft

#endif // __RR_FACTORIZATION_H__
