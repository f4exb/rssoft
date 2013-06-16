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
#include "RR_Factorization.h"
#include "GFq.h"
#include "GFq_Polynomial.h"

namespace rssoft
{

// ================================================================================================
RR_Node::RR_Node(RR_Node *_parent,
		const gf::GFq_BivariatePolynomial& _Q,
		unsigned int _id) :
	parent(_parent),
	Q(_Q),
	id(_id)
{
	if (_parent == 0)
	{
		degree = -1;
	}
	else
	{
		degree = _parent->get_degree()+1;
	}
}

// ================================================================================================
RR_Factorization::RR_Factorization(const gf::GFq& _gf, unsigned int _k) :
		gf(_gf),
		k(_k),
		t(0)
{

}

// ================================================================================================
RR_Factorization::~RR_Factorization()
{

}

// ================================================================================================
void RR_Factorization::init()
{
	t = 0;
	F.clear();
}

// ================================================================================================
void RR_Factorization::run(const gf::GFq_BivariatePolynomial& polynomial)
{
	RR_Node u(0, polynomial);
}

// ================================================================================================
gf::GFq_Polynomial RR_Factorization::node_run(RR_Node& rr_node)
{
	gf::GFq_Polynomial Qy = rr_node.getQ();
	std::vector<rssoft::gf::GFq_Element>& roots_y;
	Qy.rootChien(roots_y);
	std::vector<rssoft::gf::GFq_Element>::const_iterator ry_it = roots_y.begin();

	for (; ry_it != roots_y.end(); ++ry_it)
	{
		if (!rr_node.is_in_ry_set(*ry_it))
		{
			rr_node.add_ry(*ry_it);
		}
	}
}

} // namespace rssoft

