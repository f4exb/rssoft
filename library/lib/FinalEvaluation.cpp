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
#include "FinalEvaluation.h"
#include "GFq_Polynomial.h"
#include "EvaluationValues.h"
#include "ReliabilityMatrix.h"
#include "RSSoft_Exception.h"
 
#include <algorithm>
#include <utility>
#include <functional>
#include <iomanip>
#include <cmath>
 
namespace rssoft
{
 
// ================================================================================================
ProbabilityCodeword::ProbabilityCodeword()
{}
 
// ================================================================================================
ProbabilityCodeword::ProbabilityCodeword(float probability, const std::vector<gf::GFq_Symbol>& codeword) :
    std::pair<float, std::vector<gf::GFq_Symbol> >(probability, codeword)
{}
 
// ================================================================================================
ProbabilityCodeword::~ProbabilityCodeword()
{}

// ================================================================================================
void ProbabilityCodeword::print_codeword(std::ostream& os) const
{
	std::vector<rssoft::gf::GFq_Symbol>::const_iterator c_it = second.begin();
    os << "[";

    for (; c_it != second.end(); ++c_it)
    {
        if (c_it != second.begin())
        {
            os << ", ";
        }

        os <<  *c_it;
    }

    os << "]";
}

// ================================================================================================
FinalEvaluation::FinalEvaluation(const gf::GFq& _gf, const EvaluationValues& _evaluation_values) :
    gf(_gf),
    evaluation_values(_evaluation_values)
{
	std::vector<gf::GFq_Element>::const_iterator s_it = evaluation_values.get_symbols().begin();
    unsigned int i_s = 0;
    
    for (; s_it != evaluation_values.get_symbols().end(); ++s_it, i_s++)
    {
        symbol_index.insert(std::make_pair(*s_it, i_s));
    }
}

// ================================================================================================
FinalEvaluation::~FinalEvaluation()
{}

// ================================================================================================
void FinalEvaluation::run(const std::vector<gf::GFq_Polynomial>& polynomials, const ReliabilityMatrix& relmat)
{
    if (polynomials.size() == 0)
    {
        throw RSSoft_Exception("Cannot evaluate empty list of polynomials");
    }
    else if (relmat.get_nb_symbols() != gf.size()+1)
    {
        throw RSSoft_Exception("Reliability matrix number of rows is incompatible with GF size");
    }
    else if (relmat.get_message_length()!= evaluation_values.get_evaluation_points().size())
    {
        throw RSSoft_Exception("Reliability matrix number of columns is incompatible with the number of evaluation points");
    }
    else
    {
        std::vector<gf::GFq_Polynomial>::const_iterator poly_it = polynomials.begin();
        static const ProbabilityCodeword tmp_pc;
        
        for (; poly_it != polynomials.end(); ++poly_it)
        {
            codewords.push_back(tmp_pc);
            messages.push_back(tmp_pc);
            float proba_score = 0.0; // We will use log scale in dB/symbol
            std::vector<gf::GFq_Element>::const_iterator evalpt_it = evaluation_values.get_evaluation_points().begin();
            unsigned int i_pt = 0;
            
            for (; evalpt_it != evaluation_values.get_evaluation_points().end(); ++evalpt_it, i_pt++)
            {
                gf::GFq_Element eval = (*poly_it)(*evalpt_it); // Evaluate polynomial at current point
                codewords.back().get_codeword().push_back(eval.poly()); // Store the corresponding symbol in the codeword
                unsigned int& i_s = symbol_index.at(eval); // Retrieve symbol index in reliability matrix row order
                proba_score += 10.0 * log10(relmat(i_s, i_pt)); // Accumulate probability (dB)
            }
            
            // Probability score in dB is divided by the number of evaluation points used. This is an attempt to get a common metric among codes of different lengths
            codewords.back().get_probability_score() = proba_score/evaluation_values.get_evaluation_points().size(); // Store the probability score in the probability score weighted codeword
            messages.back().get_probability_score() = proba_score/evaluation_values.get_evaluation_points().size(); // Store the message with probability score.
            poly_it->get_poly_symbols(messages.back().get_codeword()); // Message is polynomial's coefficients
        }
    }
    
    // Sort codewords and messages vectors according to decreasing codewords probability score 
    std::sort(codewords.begin(), codewords.end(), std::greater<ProbabilityCodeword>());
    std::sort(messages.begin(), messages.end(), std::greater<ProbabilityCodeword>());
}

// ================================================================================================
void FinalEvaluation::print_codewords(std::ostream& os, const std::vector<ProbabilityCodeword>& words) const
{
	std::vector<ProbabilityCodeword>::const_iterator w_it = words.begin();
	unsigned int i_w = 0;

	for (; w_it != words.end(); ++w_it, i_w++)
	{
		std::streamsize prec = os.precision();
		os << "#" << i_w << ": (" << std::setprecision(1) << w_it->get_probability_score() << " dB/symbol) ";
		os.precision(prec);
		w_it->print_codeword(os);
		os << std::endl;
	}
}

} // namespace rssoft
 
