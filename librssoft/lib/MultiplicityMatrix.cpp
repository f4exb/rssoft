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

 Multiplicity Matrix class

 */
#include "MultiplicityMatrix.h"
#include "ReliabilityMatrix.h"
#include <iomanip>
#include <cmath>
 
namespace rssoft
{ 

// ================================================================================================
MultiplicityMatrix::MultiplicityMatrix(const ReliabilityMatrix& relmat, unsigned int multiplicity, bool soft_decision) :
    _nb_symbols_log2(relmat.get_nb_symbols_log2()),
    _nb_symbols(relmat.get_nb_symbols()),
    _message_length(relmat.get_message_length()),
    _cost(0)
{
    if (soft_decision)
    {
        ReliabilityMatrix w_relmat(relmat);
        unsigned int star_row, star_col;
        
        for (unsigned int s = multiplicity; s > 0; s--)
        {
            float p_star = w_relmat.find_max(star_row, star_col);
            iterator m_it = (*this)(star_row, star_col);
            
            if (m_it == end())
            {
                w_relmat(star_row, star_col) = p_star / 2;
                insert(std::make_pair(std::make_pair(star_row, star_col), 1));
                _cost += 1;
            }
            else
            {
                w_relmat(star_row, star_col) = p_star / (m_it->second+2);
                m_it->second += 1;
                _cost += m_it->second;
            }
        }
    }
    else // build for hard decision
    {
        for (unsigned int ic = 0; ic < _message_length; ic++)
        {
            float max_p = 0.0;
            unsigned int max_ir = 0;
            
            for (unsigned int ir = 0; ir < _nb_symbols; ir++)
            {
                if (relmat(ir, ic) > max_p)
                {
                    max_p = relmat(ir, ic);
                    max_ir = ir;
                }
            }
            
            insert(std::make_pair(std::make_pair(max_ir, ic), multiplicity));
        }
    }
}

// ================================================================================================
MultiplicityMatrix::MultiplicityMatrix(const ReliabilityMatrix& relmat, float lambda) :
    _nb_symbols_log2(relmat.get_nb_symbols_log2()),
    _nb_symbols(relmat.get_nb_symbols()),
    _message_length(relmat.get_message_length()),
    _cost(0)
 {
    for (unsigned int ic = 0; ic < _message_length; ic++)
    {
        for (unsigned int ir = 0; ir < _nb_symbols; ir++)
        {
            float p = floor(relmat(ir, ic) * lambda);
            
            if (p > 0.0)
            {
                unsigned int p_int = (unsigned int)(p);
                
                if (p_int > 0)
                {
                    insert(std::make_pair(std::make_pair(ir, ic), p_int));
                }
                
                _cost += p_int * (p_int + 1);
            }
        }
    }
    
    _cost /= 2;
 }
 
// ================================================================================================
MultiplicityMatrix::~MultiplicityMatrix()
{}

// ================================================================================================
unsigned int MultiplicityMatrix::operator()(unsigned int iX, unsigned int iY) const
{
    MultiplicityMatrix::const_iterator elt_it = this->find(std::make_pair(iX,iY));

    if (elt_it == end())
    {
    	return 0;
    }
    else
    {
    	return elt_it->second;
    }
}

// ================================================================================================
std::ostream& operator <<(std::ostream& os, const MultiplicityMatrix& matrix)
{
	unsigned int nb_rows = matrix.get_nb_symbols();
	unsigned int nb_cols = matrix.get_message_length();

	for (unsigned int ir=0; ir<nb_rows; ir++)
	{
		for (unsigned int ic=0; ic<nb_cols; ic++)
		{
			if (ic > 0)
			{
				os << " ";
			}
            
            MultiplicityMatrix::const_iterator elt_it = matrix.find(std::make_pair(ir,ic));
            
            if (elt_it == matrix.end())
            {
                os << std::setw(3) << 0;
            }
            else
            {
                os << std::setw(3) << elt_it->second;
            }
		}

		os << std::endl;
	}
    
    return os;
}

} // namespace rssoft

