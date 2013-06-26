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

 Reliability Matrix class

 */

#include "ReliabilityMatrix.h"
#include <iomanip>
#include <cstring>

namespace rssoft
{

// ================================================================================================
ReliabilityMatrix::ReliabilityMatrix(unsigned int nb_symbols_log2, unsigned int message_length) :
		_nb_symbols_log2(nb_symbols_log2),
		_nb_symbols(1<<nb_symbols_log2),
		_message_length(message_length),
		_message_symbol_count(0)
{
	_matrix = new float[_nb_symbols*_message_length];

	for (unsigned int i=0; i<_nb_symbols*_message_length; i++)
	{
		_matrix[i] = 0.0;
	}
}

// ================================================================================================
ReliabilityMatrix::ReliabilityMatrix(const ReliabilityMatrix& relmat) : 
		_nb_symbols_log2(relmat.get_nb_symbols_log2()),
		_nb_symbols(relmat.get_nb_symbols()),
		_message_length(relmat.get_message_length()),
		_message_symbol_count(0)
{
    _matrix = new float[_nb_symbols*_message_length];
    memcpy((void *) _matrix, (void *) relmat.get_raw_matrix(), _nb_symbols*_message_length*sizeof(float));
}

// ================================================================================================
ReliabilityMatrix::~ReliabilityMatrix()
{
	delete[] _matrix;
}

// ================================================================================================
void ReliabilityMatrix::enter_symbol_data(float *symbol_data)
{
	if (_message_symbol_count < _message_length)
	{
		memcpy((void *) &_matrix[_message_symbol_count*_nb_symbols], (void *) symbol_data, _nb_symbols*sizeof(float));
		_message_symbol_count++;
	}
}

// ================================================================================================
void ReliabilityMatrix::enter_symbol_data(unsigned int message_symbol_index, float *symbol_data)
{
	if (message_symbol_index < _message_length)
	{
		memcpy((void *) &_matrix[message_symbol_index*_nb_symbols], (void *) symbol_data, _nb_symbols*sizeof(float));
	}
}

// ================================================================================================
void ReliabilityMatrix::enter_erasure()
{
	if (_message_symbol_count < _message_length)
	{
        for (unsigned int i=0; i<_nb_symbols; i++)
        {
            _matrix[_message_symbol_count*_nb_symbols + i] = 0.0;
        }
        
        _message_symbol_count++;
    }    
}

// ================================================================================================
void ReliabilityMatrix::enter_erasure(unsigned int message_symbol_index)
{
	if (message_symbol_index < _message_length)
	{
        for (unsigned int i=0; i<_nb_symbols; i++)
        {
            _matrix[message_symbol_index*_nb_symbols + i] = 0.0;
        }
    }
}

// ================================================================================================
void ReliabilityMatrix::normalize()
{
	float col_sum = 0;
	float last_col_sum;

	for (unsigned int ic = 0; ic < _message_length+1; ic++)
	{
		last_col_sum = col_sum;
		col_sum = 0;

		for (unsigned int ir = 0; ir < _nb_symbols; ir++)
		{
			if (ic < _message_length)
			{
				col_sum += _matrix[ic*_nb_symbols + ir];
			}

			if (ic > 0)
			{
                if (last_col_sum != 0.0)
                {
                    _matrix[(ic-1)*_nb_symbols + ir] /= last_col_sum;
                }
			}
		}
	}
}

// ================================================================================================
float ReliabilityMatrix::find_max(unsigned int& i_row, unsigned int& i_col) const
{
    float max = 0.0;

    for (unsigned int ic = 0; ic < _message_length; ic++)
    {
        for (unsigned int ir = 0; ir < _nb_symbols; ir++)
        {
            if (_matrix[ic*_nb_symbols + ir] > max)
            {
                max = _matrix[ic*_nb_symbols + ir];
                i_row = ir;
                i_col = ic;
            }
        }
    }
    
    return max;
}

// ================================================================================================
std::ostream& operator <<(std::ostream& os, const ReliabilityMatrix& matrix)
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
			os << std::fixed << std::setw(8) << std::setprecision(6) << matrix(ir, ic);
		}

		os << std::endl;
	}
    
    return os;
}

}
