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

 Final evaluation of result polynomials to get possible codewords

 */
#ifndef __FINAL_EVALUATION_H__
#define __FINAL_EVALUATION_H__

#include "GFq.h"
#include "GFq_Element.h"
#include "GF_Utils.h"
#include <vector>
#include <map>

namespace rssoft
{
namespace gf
{
class GFq_Polynomial;
}

class EvaluationValues;
class ReliabilityMatrix;

/**
 * \brief Probability score weighted codeword 
 */
class ProbabilityCodeword : public std::pair<float, std::vector<gf::GFq_Symbol> >
{
public:
    /**
     * Constructs an empty object. Probability score is given later and codeword is constructed later
     * \param probability Probability score of this codeword
     */
    ProbabilityCodeword();
    
    /**
     * Constructs with the codeword as the copy of the given codeword
     */
    ProbabilityCodeword(float probability, const std::vector<gf::GFq_Symbol>& codeword);
    
    /**
     * Destructor. Nothing special
     */
    ~ProbabilityCodeword();
    
    /**
     * Add a new symbol to the codeword
     */
    void add_symbol(const gf::GFq_Symbol& symbol)
    {
        second.push_back(symbol);
    }
    
    /**
     * Get probability score R/O
     */
    float get_probability_score() const
    {
        return first;
    }
    
    /**
     * Get probability score R/W
     */
    float& get_probability_score()
    {
        return first;
    }
    
    /** 
     * Get the codeword
     */
    const std::vector<gf::GFq_Symbol>& get_codeword() const
    {
        return second;
    }
    
    /** 
     * Get the codeword R/W
     */
    std::vector<gf::GFq_Symbol>& get_codeword()
    {
        return second;
    }

    /** 
     * Order codewords according to their increasing probability score
     */
    bool operator<(const ProbabilityCodeword& other) const
    {
        return first < other.first;
    }
    
    /** 
     * Order codewords according to their decreasing probability score
     */
    bool operator>(const ProbabilityCodeword& other) const
    {
        return first > other.first;
    }
    
    /**
     * Print codeword to an output stream
     */
    void print_codeword(std::ostream& os) const;
};

/**
 * Print a probability score weighted codeword to an output stream
 */
std::ostream& operator <<(std::ostream& os,	const ProbabilityCodeword& scored_codeword);


class FinalEvaluation
{
public:
    /**
     * Constructor
     * \param _gf Galois Field in use
     * \param _evaluation_values Evaluation X,Y values used for coding
     */
    FinalEvaluation(const gf::GFq& _gf, const EvaluationValues& _evaluation_values);
        
    /**
     * Destructor. Nothing special
     */
    ~FinalEvaluation();
    
    /**
     * Runs one evaluation for the given polynomials
     */
    void run(const std::vector<gf::GFq_Polynomial>& polynomials, const ReliabilityMatrix& relmat);

    /**
     * Get the best probability scoring codeword
     */
    const std::vector<gf::GFq_Symbol>& get_best_codeword() const
    {
        return codewords.begin()->get_codeword();
    }
     
    /**
     * Get all codewords with their probability score
     */
    const std::vector<ProbabilityCodeword>& get_codewords() const
    {
        return codewords;
    }

    /**
     * Get the best probability scoring message
     */
    const std::vector<gf::GFq_Symbol>& get_best_message() const
    {
        return messages.begin()->get_codeword();
    }
     
    /**
     * Get all messages with their probability score
     */
    const std::vector<ProbabilityCodeword>& get_messages() const
    {
        return messages;
    }

    /**
     * Print a codeword or message word to an output stream
     */
    void print_codewords(std::ostream& os, const std::vector<ProbabilityCodeword>& words) const;

protected:
    const gf::GFq& gf; //!< Galois Field in use
    const EvaluationValues& evaluation_values; //!< Evaluation X,Y values used for coding
    std::map<gf::GFq_Element, unsigned int> symbol_index; //!< Symbol index in reliability matrix row order
    std::vector<ProbabilityCodeword> codewords; //!< The codewords (overriden at each run)
    std::vector<ProbabilityCodeword> messages; //!< The encoded messages (overriden at each run)
};

} // namespace rssoft

#endif // __FINAL_EVALUATION_H__
 
 
