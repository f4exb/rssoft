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

 Reliability Matrix class. 
 Analog data is entered first then the normalization method is called to get the actual reliability data (probabilities).

 */

#ifndef __RELIABILITY_MATRIX_H__
#define __RELIABILITY_MATRIX_H__

#include <iostream>

namespace rssoft
{

/**
 * \brief Reliability Matrix class. Analog data is entered first then the normalization method is called to get the actual reliability data (probabilities).
 */
class ReliabilityMatrix
{
public:
	/**
	 * Constructor
	 * \param nb_symbols_log2 Log2 of the number of symbols used (number of symbols is a power of two)
	 * \param message_length Length of one message block to be decoded
	 */
	ReliabilityMatrix(unsigned int nb_symbols_log2, unsigned int message_length);
    
    /**
     * Copy Constructor
     */
    ReliabilityMatrix(const ReliabilityMatrix& relmat);

	/**
	 * Destructor. Frees the matrix storage.
	 */
	~ReliabilityMatrix();

	/**
	 * Enter one more symbol position data
	 * \param symbol_data Pointer to symbol data array. There must be nb_symbol values corresponding to the relative reliability of each symbol for the current symbol position in the message
	 */
	void enter_symbol_data(float *symbol_data);

	/**
	 * Enter symbol position data at given message symbol position
	 * \param message_symbol_index Position of the symbol in the message
	 * \param symbol_data Pointer to symbol data array. There must be nb_symbol values corresponding to the relative reliability of each symbol for the current symbol position in the message
	 */
	void enter_symbol_data(unsigned int message_symbol_index, float *symbol_data);
    
    /**
     * Enter an erasure at current symbol position. This is done by zeroing out the corresponding column in the matrix thus neutralizing it for further multiplicity calculation.
     */
    void enter_erasure();

    /**
     * Enter an erasure at a given symbol position. This is done by zeroing out the corresponding column in the matrix thus neutralizing it for further multiplicity calculation.
     */
    void enter_erasure(unsigned int message_symbol_index);

	/**
	 * Normalize each column so that values represent an a posteriori probability i.e. sum of each column is 1.0
	 */
	void normalize();

	/**
	 * Resets the message symbol counter
	 */
	void reset_message_symbol_count()
	{
		_message_symbol_count = 0;
	}

	/**
	 * Get the log2 of the number of symbols (i.e. rows)
	 */
	unsigned int get_nb_symbols_log2() const
	{
		return _nb_symbols_log2;
	}

	/**
	 * Get the number of symbols (i.e. rows)
	 */
	unsigned int get_nb_symbols() const
	{
		return _nb_symbols;
	}

	/**
	 * Get the number of message symbols (i.e. columns)
	 */
	unsigned int get_message_length() const
	{
		return _message_length;
	}

	/**
	 * Operator to get the value at row i column j. Read-write version.
	 */
	float& operator()(unsigned int i_row, unsigned int i_col)
	{
		return _matrix[_nb_symbols*i_col + i_row];
	}
    
	/**
	 * Operator to get the value at row i column j. Read-only version.
	 */
	const float& operator()(unsigned int i_row, unsigned int i_col) const
	{
		return _matrix[_nb_symbols*i_col + i_row];
	}
    
    /**
     * Get a pointer to matrix storage
     */
    const float *get_raw_matrix() const
    {
        return _matrix;
    }
    
    /**
     * Finds the maximum value in the matrix
     */
    float find_max(unsigned int& i_row, unsigned int& i_col) const;

	/**
	 * Prints a reliability matrix to an output stream
	 */
	friend std::ostream& operator <<(std::ostream& os, const ReliabilityMatrix& matrix);
    

protected:
	unsigned int _nb_symbols_log2;
	unsigned int _nb_symbols;
	unsigned int _message_length;
	unsigned int _message_symbol_count; //!< incremented each time a new message symbol data is entered
	float *_matrix; //!< The reliability matrix stored column first
};


}

#endif // __RELIABILITY_MATRIX_H__
