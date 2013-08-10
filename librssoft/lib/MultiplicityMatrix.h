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
 #ifndef __MULTIPLICITY_MATRIX_H__
 #define __MULTIPLICITY_MATRIX_H__
 
 #include <utility>
 #include <map>
 #include <iostream>

namespace rssoft
{

class ReliabilityMatrix;

/**
 * \brief Ordering of elements in the sparse matrix according to the column first order. Indexes are pairs of (row, column) indexes
 */
class MultiplicityMatrix_SparseOrdering
{
public:
    MultiplicityMatrix_SparseOrdering() {}
    ~MultiplicityMatrix_SparseOrdering() {}
    
    bool operator()(const std::pair<unsigned int, unsigned int>& index1, const std::pair<unsigned int, unsigned int>& index2) const
    {
        if (index1.second == index2.second) // same column
        {
            return index1.first < index2.first; // row order
        }
        else
        {
            return index1.second < index2.second; // column order
        }
    }
};

/**
 * \brief Multiplicity matrix corresponding to a reliability matrix. It is implemented as a map representing a sparse matrix where elements are indexed by 
 * their (row, column) numbers pair. Once constructed it is normally only used to be traversed column first during the Interpolation algorithm.
 * The ordering of indexes reflects the column first order.
 */
class MultiplicityMatrix : public std::map<std::pair<unsigned int, unsigned int>, unsigned int, MultiplicityMatrix_SparseOrdering>
{
public:
    /**
     * Iterator used to traverse matrix for read only operations. Has explicit methods for indexes and value.
     */
    class traversing_iterator : public std::map<std::pair<unsigned int, unsigned int>, unsigned int, MultiplicityMatrix_SparseOrdering>::const_iterator
    {
    public:
    	/**
    	 * Constructs iterator from the parent map const iterator. Thus it can be initialized with something like:
    	 * traversing_iterator it(matrix.begin());
    	 */
    	traversing_iterator(std::map<std::pair<unsigned int, unsigned int>, unsigned int, MultiplicityMatrix_SparseOrdering>::const_iterator const &c) :
    		std::map<std::pair<unsigned int, unsigned int>, unsigned int, MultiplicityMatrix_SparseOrdering>::const_iterator(c) {}

    	/**
    	 * Return the index in X which is the column coordinate
    	 */
    	unsigned int iX()
    	{
    		std::map<std::pair<unsigned int, unsigned int>, unsigned int, MultiplicityMatrix_SparseOrdering>::const_iterator *c_it = this;
    		return (*c_it)->first.second;
    	}

    	/**
    	 * Return the index in Y which is the row coordinate
    	 */
    	unsigned int iY()
    	{
    		std::map<std::pair<unsigned int, unsigned int>, unsigned int, MultiplicityMatrix_SparseOrdering>::const_iterator *c_it = this;
    		return (*c_it)->first.first;
    	}

    	/**
    	 * Return the multiplicity value
    	 */
    	unsigned int multiplicity()
    	{
    		std::map<std::pair<unsigned int, unsigned int>, unsigned int, MultiplicityMatrix_SparseOrdering>::const_iterator *c_it = this;
    		return (*c_it)->second;
    	}
    };

    /**
     * Constructs a new multiplicity matrix. Uses long construction algorithm for soft decision.
     * \param relmat Reliability matrix to build the multiplicity matrix from
     * \param multiplicity For soft decision: target global multiplicity of interpolation points. For hard decision: multiplicity at each point
     * \param soft_decision True to build the matrix for soft decision decoding (default) 
     *                      else to build the matrix for hard decision list decoding with specified multiplicity at each point
     */
    MultiplicityMatrix(const ReliabilityMatrix& relmat, unsigned int multiplicity, bool soft_decision=true);

    /**
     * Constructs a new multiplicity matrix. Uses short construction algorithm.
     * \param relmat Reliability matrix to build the multiplicity matrix from
     * \param lambda Multiplicative constant
     */
    MultiplicityMatrix(const ReliabilityMatrix& relmat, float lambda);
    
    /**
     * Destructor
     */
    ~MultiplicityMatrix();
    
    /**
     * Get multiplicity matrix cost
     * \return Multiplicity matrix cost
     */
    unsigned int cost() const
    {
        return _cost;
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
	 * Operator to get value at row i column j.
	 */
	unsigned int operator()(unsigned int i_row, unsigned int i_col) const;

	/**
	 * Prints a multiplicity matrix to an output stream
	 */
	friend std::ostream& operator <<(std::ostream& os, const MultiplicityMatrix& matrix);
    

protected:
    /**
     * Multiplicity matrix const iterator. The multiplicity matrix is a map of multiplicities (unsigned int) indexed by (row, column) pair
     * Thus if it is an iterator of MultiplicityMatrix::const_iterator type:
     *   - it->first.first is the row index or symbol in the alphabet
     *   - it->first.second is the column index or symbol index in the message and point of evaluation in GFq
     *   - it->second is the multiplicity value
     */
    typedef std::map<std::pair<unsigned int, unsigned int>, unsigned int, MultiplicityMatrix_SparseOrdering>::const_iterator const_iterator;

    /**
     * Multiplicity matrix iterator. The multiplicity matrix is a map of multiplicities (unsigned int) indexed by (row, column) pair
     * Thus if it is an iterator of MultiplicityMatrix::const_iterator type:
     *   - it->first.first is the row index or symbol in the alphabet
     *   - it->first.second is the column index or symbol index in the message and point of evaluation in GFq
     *   - it->second is the multiplicity value
     */
    typedef std::map<std::pair<unsigned int, unsigned int>, unsigned int, MultiplicityMatrix_SparseOrdering>::iterator iterator;

	/**
	 * Operator to get iterator at row i column j. Read-write version.
	 */
	iterator operator()(unsigned int i_row, unsigned int i_col)
    {
        return find(std::make_pair(i_row, i_col));
    }

	unsigned int _nb_symbols_log2; //!< log2 of the number of symbols in the alphabet
	unsigned int _nb_symbols; //!< number of symbols in the alphabet
	unsigned int _message_length; //!< message or block length
    unsigned int _cost; //!< Multiplicity matrix cost
};
 
} // namespace rssoft

 #endif // __MULTIPLICITY_MATRIX_H__
 
