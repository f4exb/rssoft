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

	 Full coding/decoding test with direct AWGN

*/

#include "GFq.h"
#include "GF2_Element.h"
#include "GF2_Polynomial.h"
#include "GF_Utils.h"
#include "EvaluationValues.h"
#include "ReliabilityMatrix.h"
#include "MultiplicityMatrix.h"
#include "GSKV_Interpolation.h"
#include "RR_Factorization.h"
#include "FinalEvaluation.h"
#include "RS_Encoding.h"
#include "URandom.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <getopt.h>
#include <boost/lexical_cast.hpp>

// ================================================================================================
// template to extract information from getopt more easily
template<typename TOpt, typename TField> bool extract_option(TField& field, char short_option)
{
    TOpt option_value;
    
    try
    {
        option_value = boost::lexical_cast<TOpt>(optarg);
        field = option_value;
        return true;
    }
    catch (boost::bad_lexical_cast &)
    {
        std::cout << "wrong argument for -" << short_option << ": " << optarg << " leave default (" << field << ")";
        std::cout << std::endl;
        return false;
    }
}

// ================================================================================================
struct Options
{
public:
    Options() :
        make_noise(false),
        snr_dB(0),
        m(3),
        k(5),
        global_multiplicity(1<<3),
        verbosity(0),
        print_seed(false),
        print_stats(false),
        seed(0),
        has_seed(false),
        print_sagemath(false),
        iterations(1),
        nb_erasures(0),
        _indicator_int(0)
    {
        // http://theory.cs.uvic.ca/gen/poly.html
        rssoft::gf::GF2_Element pp_gf8[4]   = {1,1,0,1};
        rssoft::gf::GF2_Element pp_gf16[5]  = {1,0,0,1,1};
        rssoft::gf::GF2_Element pp_gf32[6]  = {1,0,0,1,0,1};
        rssoft::gf::GF2_Element pp_gf64[7]  = {1,0,0,0,0,1,1};
        rssoft::gf::GF2_Element pp_gf128[8] = {1,0,0,0,0,0,1,1};
        rssoft::gf::GF2_Element pp_gf256[9] = {1,0,0,0,1,1,1,0,1};
    
        ppolys.push_back(rssoft::gf::GF2_Polynomial(4,pp_gf8));
        ppolys.push_back(rssoft::gf::GF2_Polynomial(5,pp_gf16));
        ppolys.push_back(rssoft::gf::GF2_Polynomial(6,pp_gf32));
        ppolys.push_back(rssoft::gf::GF2_Polynomial(7,pp_gf64));
        ppolys.push_back(rssoft::gf::GF2_Polynomial(8,pp_gf128));
        ppolys.push_back(rssoft::gf::GF2_Polynomial(9,pp_gf256));
    }

    const rssoft::gf::GF2_Polynomial& get_ppoly() const
    {
        return ppolys[m-3];
    }
    
    bool get_options(int argc, char *argv[]);
    
    bool make_noise;
    float snr_dB;
    unsigned int m;
    unsigned int k;
    unsigned int global_multiplicity;
    unsigned int verbosity;
    bool print_seed;
    bool print_stats;
    unsigned int seed;
    bool has_seed;
    bool print_sagemath; //!< Print input for Sage Math script
    unsigned int iterations; //!< Maximum number of retry iterations
    unsigned int nb_erasures; //!< Number of erasures
    int _indicator_int;
private:
    std::vector<rssoft::gf::GF2_Polynomial> ppolys;
};

// ================================================================================================
struct StatOutput
{
public:
	StatOutput() :
		codeword_average_score(0.0),
		snr_dB(0.0),
		nb_hard_errors(0),
		nb_erasures(0),
		nb_results_when_found(0),
		result_order_when_found(0),
		found(false),
		nb_iterations(0),
		nb_false_results(0),
		max_multiplicity(0),
		max_matrix_cost(0)
	{}

	friend std::ostream& operator<<(std::ostream& os, const StatOutput& polynomial);

	float codeword_average_score;
	float snr_dB;
	unsigned int nb_hard_errors;
    unsigned int nb_erasures;
	unsigned int nb_results_when_found;
	unsigned int result_order_when_found;
	bool found;
	unsigned int nb_iterations;
	unsigned int nb_false_results;
	unsigned int max_multiplicity;
	unsigned int max_matrix_cost;
};


// ================================================================================================
std::ostream& operator<<(std::ostream& os, const StatOutput& st)
{
	os << st.snr_dB << ","
		<< st.codeword_average_score << ","
		<< st.nb_hard_errors << ","
		<< st.nb_erasures << ","
		<< st.found << ","
		<< st.nb_results_when_found << ","
		<< st.result_order_when_found << ","
		<< st.nb_iterations << ","
		<< st.nb_false_results << ","
		<< st.max_multiplicity << ","
		<< st.max_matrix_cost;
	return os;
}


// ================================================================================================
bool Options::get_options(int argc, char *argv[])
{
    int c;
    bool status = true;

    while (true)
    {
        static struct option long_options[] =
        {
            // these options set a flag
            {"print-seed", no_argument, &_indicator_int, 1},
            {"print-stats", no_argument, &_indicator_int, 1},
            {"sagemath", no_argument, &_indicator_int, 1},
            // these options do not set a flag
            {"snr", required_argument, 0, 'n'},        
            {"log2-n", required_argument, 0, 'm'},      
            {"k", required_argument, 0, 'k'},              
            {"global-multiplicity", required_argument, 0, 'M'},
            {"verbosity", required_argument, 0, 'v'},              
            {"seed", required_argument, 0, 's'},              
            {"nb-iterations-max", required_argument, 0, 'i'},
            {"nb-erasures", required_argument, 0, 'e'},
        };    
        
        int option_index = 0;
        c = getopt_long (argc, argv, "n:m:k:M:v:s:i:e:", long_options, &option_index);
        
        if (c == -1) // end of options
        {
            break;
        }

        switch(c)
        {
            case 0: // set flag
                if (strcmp("print-seed", long_options[option_index].name) == 0)
                {
                    print_seed = true;
                }
                if (strcmp("print-stats", long_options[option_index].name) == 0)
                {
                    print_stats = true;
                }
                if (strcmp("sagemath", long_options[option_index].name) == 0)
                {
                    print_sagemath = true;
                }
                _indicator_int = 0;
                break;
            case 'n':
                make_noise = true;
                status = extract_option<double, float>(snr_dB, 'n');
                break;
            case 'm':
                status = extract_option<int, unsigned int>(m, 'm');
                break;
            case 'k':
                status = extract_option<int, unsigned int>(k, 'k');
                break;
            case 'M':
                status = extract_option<int, unsigned int>(global_multiplicity, 'M');
                break;
            case 'v':
                status = extract_option<int, unsigned int>(verbosity, 'v');
                break;
            case 's':
                status = extract_option<int, unsigned int>(seed, 's');
                has_seed = true;
                break;
            case 'i':
                status = extract_option<int, unsigned int>(iterations, 'i');
                break;
            case 'e':
                status = extract_option<int, unsigned int>(nb_erasures, 'e');
                break;
            case '?':
                status = false;
                break;
        }    
        
    }
    
    if (status)
    {
        unsigned int n = (1<<m) - 1;
        
        if ((m < 3) || (m > 8))
        {
            std::cout << "Not implemented for GF(2^" << m << ") fields" << std::endl;
            status = false;
        }
    
        if (k > (n-2))
        {
            std::cout << "Cannot work with RS(" << n << "," << k << ")" << std::endl;
            status = false;
        }
        else if (k < 2)
        {
            std::cout << "Cannot work with RS(" << n << "," << k << ")" << std::endl;
            status = false;
        }

        if (global_multiplicity < n)
        {
        	std::cout << "Global multiplicity must be at least " << n << std::endl;
        	status = false;
        }
        
        if (nb_erasures > n-2)
        {
            std::cout << "The number of erasures (" << nb_erasures << ") cannot exceed the number of symbols - 2 (" << n-2 << ")" << std::endl;
        }
    }
    
    return status;
}


// ================================================================================================
int main(int argc, char *argv[])
{
    Options options;
    
    if (options.get_options(argc, argv))
    {
        unsigned int q = (1<<options.m);
        unsigned int n = q - 1;
        double std_dev  = 1.0 / pow(10.0, (options.snr_dB/10.0)); // Standard deviation for power AWGN
        StatOutput stat_output;
        std::set<unsigned int> erased_indexes;

        rssoft::gf::GFq gfq(options.m, options.get_ppoly());
        
        URandom ur;
        
        if (options.has_seed)
        {
            ur.set_seed(options.seed);
        }
        
        if (options.print_seed)
        {
            unsigned int seed = ur.rand_uword();
            std::cout << "Seed = " << seed << std::endl;
            ur.set_seed(seed);
        }
        
        std::vector<rssoft::gf::GFq_Symbol> message;
        
        for (unsigned int i=0; i<options.k; i++)
        {
            message.push_back(ur.rand_int(q));
        }
        
        std::cout << "Message : (k=" << message.size() << ") ";
        rssoft::gf::print_symbols_vector(std::cout, message);
        std::cout << std::endl;
        rssoft::EvaluationValues evaluation_values(gfq); // use default
        rssoft::RS_Encoding rs_encoding(gfq, options.k, evaluation_values);
        std::vector<rssoft::gf::GFq_Symbol> codeword;
        std::vector<unsigned int> row_indexes;

        rs_encoding.run(message, codeword);

        std::cout << "Codeword: (n=" << codeword.size() << ") ";
        rssoft::gf::print_symbols_vector(std::cout, codeword);
        std::cout << std::endl;

        // Prepare set of erased symbol indexes:
        
        std::set<unsigned int>::iterator s_it;
        
        while (erased_indexes.size() < options.nb_erasures)
        {
            unsigned int i_s = ur.rand_int(n);
            erased_indexes.insert(i_s);
        }
        
        std::cout << "Erasures: (n=" << codeword.size() << ") ";
        rssoft::gf::print_symbols_and_erasures(std::cout, codeword, erased_indexes);
        std::cout << std::endl;

        // Simulate reception behind noisy channel. Reliability matrix is created with noisy power samples.

        rssoft::ReliabilityMatrix mat_Pi(options.m,n);
        std::vector<rssoft::gf::GFq_Symbol> hard_decision;
        unsigned int hard_decision_errors = 0;
        float *mat_Pi_col = new float[q];

        for (unsigned int c=0; c<n; c++)
        {
        	float max_pwr = 0;
        	unsigned int r_max = 0;

            if (erased_indexes.find(c) == erased_indexes.end()) // symbol not erased
            {
                for (unsigned int r=0; r<q; r++)
                {
                    if (evaluation_values.get_y_values()[r] == codeword[c]) // evaluation point
                    {
                        row_indexes.push_back(r);
                        mat_Pi_col[r] = 1.0 + (options.make_noise ? std_dev * ur.rand_gaussian() : 0.0);
                        mat_Pi_col[r] *= mat_Pi_col[r];
                    }
                    else
                    {
                        mat_Pi_col[r] = 0.0 + (options.make_noise ? std_dev * ur.rand_gaussian() : 0.0);
                        mat_Pi_col[r] *= mat_Pi_col[r];
                    }

                    if (mat_Pi_col[r] > max_pwr)
                    {
                        max_pwr = mat_Pi_col[r];
                        r_max = r;
                    }
                }

                mat_Pi.enter_symbol_data(mat_Pi_col);
                hard_decision.push_back(evaluation_values.get_y_values()[r_max].poly());

                if (hard_decision.back() != codeword[c])
                {
                    hard_decision_errors++;
                }
            }
            else // erased symbol
            {
                mat_Pi.enter_erasure();
                row_indexes.push_back(0);
                hard_decision.push_back(0);
            }
        }

        delete[] mat_Pi_col;

        std::cout << "Hard-dec: (n=" << codeword.size() << ") ";
        rssoft::gf::print_symbols_and_erasures(std::cout, hard_decision, erased_indexes);
        std::cout << std::endl;
        std::cout << " -> " << hard_decision_errors << " errors, " << options.nb_erasures << " erasures: " << ((2*hard_decision_errors)+options.nb_erasures < (n-options.k) ? "correctable" : "uncorrectable") << " with hard decision" << std::endl;

        if (options.verbosity > 0)
        {
        	std::cout << "Power matrix:" << std::endl;
        	std::cout << mat_Pi;
        	std::cout << std::endl;
        }

        mat_Pi.normalize();
        float codeword_score = 0.0;
        unsigned int codeword_count = 0;
        float best_score, worst_score;
        bool first_round = true;

        for (unsigned int c=0; c<n; c++)
        {
            if (erased_indexes.find(c) == erased_indexes.end()) // symbol not erased
            {
                float score = 10.0 * log10(mat_Pi(row_indexes[c],c));

                if (first_round)
                {
                    best_score = score;
                    worst_score = score;
                    first_round = false;
                }
                else
                {
                    if (score > best_score)
                    {
                        best_score = score;
                    }
                    if (score < worst_score)
                    {
                        worst_score = score;
                    }
                }

                codeword_score += score;
                codeword_count++;
            }
        }

        stat_output.snr_dB = options.snr_dB;
        stat_output.codeword_average_score = codeword_score / codeword_count;
        stat_output.nb_hard_errors = hard_decision_errors;
        stat_output.nb_erasures = options.nb_erasures;

        std::cout << "Codeword score: " << codeword_score / codeword_count << " dB/symbol (best = " << best_score << ", worst = " << worst_score << ")" << std::endl;
        bool found = false;
        unsigned int global_multiplicity = options.global_multiplicity;

        for (unsigned int ni=1; (ni<=options.iterations) && (!found); ni++)
        {
   			std::cout << std::endl;
			rssoft::MultiplicityMatrix mat_M(mat_Pi, global_multiplicity);

			if (options.verbosity > 0)
			{
				std::cout << "Multiplicity matrix:" << std::endl;
				std::cout << mat_M;
				std::cout << std::endl;
			}

			unsigned int mm_cost = mat_M.cost();
			std::cout << "Multiplicity matrix cost is " << mm_cost << std::endl;

			rssoft::GSKV_Interpolation gskv(gfq, options.k, evaluation_values);
			rssoft::RR_Factorization rr(gfq, options.k);
			gskv.set_verbosity(options.verbosity);
			rr.set_verbosity(options.verbosity);

			const rssoft::gf::GFq_BivariatePolynomial& Q = gskv.run(mat_M);
			std::cout << "Q(X,Y) = " << Q << std::endl;

			if (Q.is_in_X())
			{
				std::cout << "Interpolation polynomial is in X only and is not factorizable. Hence no solutions" << std::endl;
			}
			else
			{
				std::vector<rssoft::gf::GFq_Polynomial>& res_polys = rr.run(Q);

				std::cout << res_polys.size() << " result(s)" << std::endl;

				if (res_polys.size() > 0)
				{
					std::vector<rssoft::gf::GFq_Polynomial>::iterator respoly_it = res_polys.begin();
					unsigned int i=0;
					for (; respoly_it != res_polys.end(); ++respoly_it, i++)
					{
						respoly_it->set_alpha_format(true);
						std::cout << "F" << i << "(X) = " << *respoly_it << std::endl;
					}

					rssoft::FinalEvaluation final_evaluation(gfq, options.k, evaluation_values);
					final_evaluation.run(res_polys, mat_Pi);
					std::cout << "Codewords:" << std::endl;
					final_evaluation.print_codewords(std::cout, final_evaluation.get_codewords());
					std::cout << "Messages:" << std::endl;
					const std::vector<rssoft::ProbabilityCodeword>& messages = final_evaluation.get_messages();
					final_evaluation.print_codewords(std::cout, messages);

					std::vector<rssoft::ProbabilityCodeword>::const_iterator ms_it = messages.begin();
					unsigned int i_m = 0;

					for (; ms_it != messages.end(); ++ms_it, i_m++)
					{
						if (rssoft::gf::compare_symbol_vectors(ms_it->get_codeword(), message))
						{
							std::cout << "#" << i_m << " found at iteration #" << ni << " !!!" << std::endl;
							stat_output.found = true;
							stat_output.nb_results_when_found = res_polys.size();
							stat_output.result_order_when_found = i_m;
							stat_output.max_matrix_cost = mm_cost;
							stat_output.max_multiplicity = global_multiplicity;
							found = true;
						}
						else
						{
							stat_output.nb_false_results++;
						}
					}
				}
			}

			if (options.print_sagemath)
			{
				std::cout << "    dY=" << gskv.get_dY() << std::endl;
				std::cout << "    Cm=" << mat_M.cost() << std::endl;
				std::cout << "    k=" << options.k << std::endl;
				std::vector<rssoft::gf::GFq_Element> x_values;
				std::vector<rssoft::gf::GFq_Element> y_values;
				std::vector<unsigned int> multiplicities;

				rssoft::MultiplicityMatrix::traversing_iterator m_it(mat_M.begin());

				for (; m_it != mat_M.end(); ++m_it)
				{
					x_values.push_back(evaluation_values.get_x_values()[m_it.iX()]);
					y_values.push_back(evaluation_values.get_y_values()[m_it.iY()]);
					multiplicities.push_back(m_it.multiplicity());
				}

				std::vector<rssoft::gf::GFq_Element>::const_iterator gfe_it = x_values.begin();
				std::cout << "    x=[";

				for (; gfe_it != x_values.end(); ++gfe_it)
				{
					if (gfe_it != x_values.begin())
					{
						std::cout << ",";
					}

					std::cout << *gfe_it;
				}

				std::cout << "]" << std::endl;
				gfe_it = y_values.begin();
				std::cout << "    y=[";

				for (; gfe_it != y_values.end(); ++gfe_it)
				{
					if (gfe_it != y_values.begin())
					{
						std::cout << ",";
					}

					std::cout << *gfe_it;
				}

				std::cout << "]" << std::endl;
				std::vector<unsigned int>::const_iterator mul_it = multiplicities.begin();
				std::cout << "    m=[";

				for (; mul_it != multiplicities.end(); ++mul_it)
				{
					if (mul_it != multiplicities.begin())
					{
						std::cout << ",";
					}

					std::cout << *mul_it;
				}

				std::cout << "]" << std::endl;
			}

			if (!found)
			{
				stat_output.max_matrix_cost = mm_cost;
				stat_output.max_multiplicity = global_multiplicity;
			}

			global_multiplicity = mm_cost;
			stat_output.nb_iterations = ni;
        } // retry iterations

        if (options.print_stats)
        {
            std::cout << std::endl;
        	std::cout << "#RES: " << stat_output << std::endl;
        }

        return 0;
    }
    else
    {
        std::cout << "Wrong options" << std::endl;
        return -1;
    }
}
