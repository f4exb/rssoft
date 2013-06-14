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

 Guruswami-Sudan-Koetter-Vardy interpolation class for soft decision
 decoding

 */
#include "GSKV_Interpolation.h"
#include "RSSoft_Exception.h"
#include "GFq.h"
#include "MultiplicityMatrix.h"
#include <cmath>

namespace rssoft
{

// ================================================================================================
GSKV_Interpolation::GSKV_Interpolation(const gf::GFq& _gf, unsigned int _k) :
		gf(_gf),
		k(_k)
{
	if (k < 2)
	{
		throw RSSoft_Exception("k parameter must be at least 2");
	}

	// default X interpolation values initialized as increasing powers of alpha starting at a^0 = 1
	// default Y interpolation values initialized as 0 followed by increasing powers of alpha starting at a^0 = 1

	y_values.push_back(gf::GFq_Element(gf, 0));

	for (unsigned int i=0; i < gf.size(); i++)
	{
		x_values.push_back(gf::GFq_Element(gf, gf.alpha(i)));
		y_values.push_back(gf::GFq_Element(gf, gf.alpha(i)));
	}
}

// ================================================================================================
GSKV_Interpolation::~GSKV_Interpolation()
{}

// ================================================================================================
void GSKV_Interpolation::run(const MultiplicityMatrix& mmat)
{
	std::pair<unsigned int, unsigned int> max_degrees = maximum_degrees(mmat);
	unsigned int dX = max_degrees.first;
	unsigned int dY = max_degrees.second;
	std::cout << "dX = " << dX << ", dY = " << dY << std::endl;

	init_G(dY);

	// outer loop on multiplicity matrix elements

	MultiplicityMatrix::traversing_iterator m_it(mmat.begin());

	for (; m_it != mmat.end(); ++m_it)
	{
		process_point(m_it.iX(), m_it.iY(), m_it.multiplicity());
	}
}

// ================================================================================================
std::pair<unsigned int, unsigned int> GSKV_Interpolation::maximum_degrees(const MultiplicityMatrix& mmat)
{
	float fdX, fdY;

	fdY = floor((1+sqrt(1+((8*mmat.cost())/(k-1))))/2.0f) - 1.0;
	fdX = floor((mmat.cost()/(fdY+1)) + ((fdY*(k-1))/2));

	return std::make_pair((unsigned int) fdX, (unsigned int) fdY);
}

// ================================================================================================
void GSKV_Interpolation::init_G(unsigned int dY)
{
	unsigned int inclod = 1;
	unsigned int lod = 0;

	for (unsigned int i=0; i<dY+1; i++)
	{
		gf::GFq_BivariatePolynomial Y_i(1, k-1);
		Y_i.init_y_pow(gf, i);
		G.push_back(Y_i);
		calcG.push_back(true);
		lodG.push_back(lod);
		inclod += k-1;
		lod += inclod;
	}
}

// ================================================================================================
void GSKV_Interpolation::process_point(unsigned int iX, unsigned int iY, unsigned int multiplicity)
{
	for (unsigned int mu = 0; mu < multiplicity; mu++)
	{
		for (unsigned int nu = 0; nu < multiplicity-mu; nu++)
		{
			process_hasse(x_values[iX], y_values[iY], mu, nu);
		}
	}
}

// ================================================================================================
void GSKV_Interpolation::process_hasse(const gf::GFq_Element& x, const gf::GFq_Element& y, unsigned int mu, unsigned int nu)
{

}

} // namespace rssoft

