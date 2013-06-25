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
#include "GFq_Element.h"
#include "GFq_BivariatePolynomial.h"
#include "EvaluationValues.h"
#include "MultiplicityMatrix.h"
#include "Debug.h"
#include <cmath>
#include <iostream>
#include <string>
#include <algorithm>

namespace rssoft
{

// ================================================================================================
GSKV_Interpolation::GSKV_Interpolation(const gf::GFq& _gf, unsigned int _k, const EvaluationValues& _evaluation_values) :
		gf(_gf),
		k(_k),
		evaluation_values(_evaluation_values),
		it_number(0),
		Cm(0),
		final_ig(0),
        verbosity(0),
        dX(0),
        dY(0),
        mcost(0)
{
	if (k < 2)
	{
		throw RSSoft_Exception("k parameter must be at least 2");
	}
}

// ================================================================================================
GSKV_Interpolation::~GSKV_Interpolation()
{}

// ================================================================================================
const gf::GFq_BivariatePolynomial& GSKV_Interpolation::run(const MultiplicityMatrix& mmat)
{
	std::pair<unsigned int, unsigned int> max_degrees = maximum_degrees(mmat);
	dX = max_degrees.first;
	dY = max_degrees.second;
	mcost = mmat.cost();

    //DebugMsg(0,verbose) << "dX = ";
    //DebugStream() << "dX = " << "toto" << std::endl;
	DEBUG_OUT(verbosity > 0, "dX = " << dX << ", dY = " << dY << std::endl);

	init_G(dY);
    it_number = 0;
    Cm = mmat.cost();

	// outer loop on multiplicity matrix elements

    DEBUG_OUT(verbosity > 0, "Loop on multiplicity matrix elements:" << std::endl);
	MultiplicityMatrix::traversing_iterator m_it(mmat.begin());

	for (; m_it != mmat.end(); ++m_it)
	{
        DEBUG_OUT(verbosity > 0, "*** New point iX = " << m_it.iX() << " iY = " << m_it.iY() << " mult = " << m_it.multiplicity() <<  std::endl);
		process_point(m_it.iX(), m_it.iY(), m_it.multiplicity());
	}

	const gf::GFq_BivariatePolynomial& Q = final_G();
	DEBUG_OUT(verbosity > 0, "Q(X,Y) = " << Q << std::endl);
    return Q;
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
		G.push_back(Y_i);
		G.back().init_y_pow(gf, i);
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
			process_hasse(evaluation_values.get_x_values()[iX], evaluation_values.get_y_values()[iY], mu, nu);
		}
	}
}

// ================================================================================================
void GSKV_Interpolation::process_hasse(const gf::GFq_Element& x, const gf::GFq_Element& y, unsigned int mu, unsigned int nu)
{
    unsigned int ig_lodmin = 0; //!< index of polynomial in G with minimal leading order
    unsigned int lodmin = 0;    //!< minimal leading order of polynomials in G
    bool first_hnn = true;
    std::vector<gf::GFq_Element> hasse_xy_G;             //!< evaluations of Hasse derivative at (x,y) for all polynomials in G
    std::vector<gf::GFq_BivariatePolynomial> G_next; //!< G list for next iteration
    std::vector<unsigned int> lodG_next;             //!< Leading orders of polynomials in G_next
    bool zero_Hasse = true;
    std::string ind("");        //!< indicator character for debug display
    
    DEBUG_OUT(verbosity > 1, "it=" << it_number << " x=" << x << " y=" << y << " mu=" << mu << " nu=" << nu << " G.size()=" << G.size() << std::endl);
    
    // Hasse derivatives calculation
    
    unsigned int ig = 0;
    std::vector<gf::GFq_BivariatePolynomial>::const_iterator it_g = G.begin();
    
    for (; it_g != G.end(); ++it_g, ig++)
    {
        if (calcG[ig]) // Polynomial is part of calculation as per Li Chen's optimization
        {
            gf::GFq_BivariatePolynomial h = dHasse(mu, nu, *it_g);
            hasse_xy_G.push_back(h(x,y));
            unsigned int wd = it_g->wdeg();
            
            if (hasse_xy_G.back().is_zero())
            {
                ind = "=";
            }
            else
            {
                zero_Hasse = false;
                ind = "!";
                
                // initialize polynomial localization variables
                if (first_hnn)
                {
                    lodmin = lodG[ig];
                    ig_lodmin = ig;
                    first_hnn = false;
                }
                
                // locate polynomial in G with minimal leading order
                if (lodG[ig] < lodmin) 
                {
                    lodmin = lodG[ig];
                    ig_lodmin = ig;
                }
            }
        }
        else // Polynomial is skipped for calculation due to Li Chen's optimization
        {
            ind = "x";
            hasse_xy_G.push_back(gf::GFq_Element(gf,0));
        }
        
        // debug print stuff
        DEBUG_OUT(verbosity > 1, ind << " G_" << it_number << "[" << ig << "] = " << *it_g << std::endl);
        
        if (calcG[ig])
        {
            DEBUG_OUT(verbosity > 1, "  D_" << it_number << "," << ig << " = " << hasse_xy_G.back() << std::endl);
            DEBUG_OUT(verbosity > 1, "  lod = " << lodG[ig] << std::endl);
        }
        else
        {
            DEBUG_OUT(verbosity > 1, "  lod = " << lodG[ig] << std::endl);
        }
    }

    DEBUG_OUT(verbosity > 1, "Minimal LOD polynomial G_" << it_number << "[" << ig_lodmin << "]" << std::endl);

    if (zero_Hasse)
    {
        DEBUG_OUT(verbosity > 1, "All Hasse derivatives are 0 so G_" << it_number+1 << " = G_" << it_number << std::endl);
    }
    else
    {
        // compute next values in G
    	it_g = G.begin();
		ig = 0;

		for (; it_g != G.end(); ++it_g, ig++)
		{
			if (calcG[ig]) // Polynomial is part of calculation as per Li Chen's optimization
			{
				if (hasse_xy_G[ig].is_zero())
				{
					G_next.push_back(*it_g); // carry over the same polynomial
					lodG_next.push_back(lodG[ig]);
				}
				else
				{
					if (ig == ig_lodmin) // Polynomial with minimal leading order
					{
						gf::GFq_BivariatePolynomial X1(1,k-1);
						X1.init_x_pow(gf,1); // X1(X,Y) = X
						G_next.push_back(hasse_xy_G[ig]*(*it_g)*(X1-x));
						unsigned int mX = it_g->lmX(); // leading monomial's X power
						unsigned int mY = it_g->lmY(); // leading monomial's Y power
						lodG_next.push_back(lodG[ig_lodmin]+(mX/(k-1))+1+mY); // new leading order by sliding one position of X powers to the right
					}
					else // other polynomials
					{
						G_next.push_back(hasse_xy_G[ig]*G[ig_lodmin]-hasse_xy_G[ig_lodmin]*(*it_g));
						lodG_next.push_back(std::max(lodG[ig],lodG[ig_lodmin]));   // new leading order is the max of the two
					}
				}

				if (lodG_next.back() > Cm)
				{
					calcG[ig] = false; // Li Chen's complexity reduction, skip polynomial processing if its lod is too big (bigger than multiplicity cost)
				}
			}
			else // Polynomial is skipped for calculation due to Li Chen's optimization
			{
				G_next.push_back(*it_g); // carry over the same polynomial
				lodG_next.push_back(lodG[ig]);
			}
		}

		// store next values if matrix cost is not reached
		if (it_number < mcost)
		G.assign(G_next.begin(), G_next.end());
		lodG.assign(lodG_next.begin(), lodG_next.end());
		it_number++;
    }
    
	DEBUG_OUT(verbosity > 1, std::endl);
}

// ================================================================================================
const gf::GFq_BivariatePolynomial& GSKV_Interpolation::final_G()
{
    bool first_g = true;
    unsigned int ig_lodmin = 0;    //!< index of polynomial in G with minimal leading order
    unsigned int lodmin = lodG[0]; //!< minimal leading order of polynomials in G

    unsigned int ig = 0;
    std::vector<gf::GFq_BivariatePolynomial>::const_iterator it_g = G.begin();

    DEBUG_OUT(verbosity > 1, "it=" << it_number << " final result" << std::endl);

    for (; it_g != G.end(); ++it_g, ig++)
    {
        /*
        if (it_g->is_in_X()) // Circumvent design flaw where an X solution only can be returned which cannot be factorized. (Is this a design flaw?)
        {
            //DEBUG_OUT(verbosity > 1, "g_" << ig << " is a polynomial in X only" << std::endl);
            std::cout << "g_" << ig << " is a polynomial in X only with LOD = " << lodG[ig] << std::endl;
            continue;
        }
        */
    
        if (first_g)
        {
            ig_lodmin = ig;
            lodmin = lodG[ig];
            first_g = false;
        }
    
    	if (lodG[ig] < lodmin)
    	{
    		lodmin = lodG[ig];
    		ig_lodmin = ig;
    	}

    	DEBUG_OUT(verbosity > 1, "o G_" << it_number << "[" << ig << "] = " << *it_g << std::endl);
    	DEBUG_OUT(verbosity > 1, "  lod = " << lodG[ig] << std::endl);
    }

    DEBUG_OUT(verbosity > 1, "Minimal LOD polynomial G_" << it_number << "[" << ig_lodmin << "]" << std::endl);
    //std::cout << "Min LOD = " << lodmin << std::endl;
    return G[ig_lodmin];
}

} // namespace rssoft

