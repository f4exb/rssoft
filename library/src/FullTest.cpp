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
        verbosity(0)
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
private:
    std::vector<rssoft::gf::GF2_Polynomial> ppolys;
};

// ================================================================================================
class URandom
{
public:
    URandom() 
    {
        rf = fopen("/dev/urandom", "r");
    }
    
    ~URandom() 
    {
        fclose(rf);
    }
    
    int rand_word() 
    {
        int ri; 
        unsigned int bytes_read = fread((char*)(&ri),sizeof(ri),1,rf);
        return ri & 0x7fffffff;
    }
    
    double rand_uniform() // GENERATE UNIFORMLY FROM [0,1)
    {
        return (double)rand_word() / (1.0+(double)0x7fffffff);
    }
    
    double rand_uniopen() // GENERATE UNIFORMLY FORM (0,1)
    {
        return (0.5+(double)rand_word()) / (1.0+(double)0x7fffffff);
    }
    
    unsigned int rand_int(unsigned int n) // GENERATE RANDOM INTEGER FROM 0, 1, ..., (n-1)
    { 
      return (unsigned int) (n * rand_uniform());
    }    
    
    double rand_gaussian() // GAUSSIAN GENERATOR.  Done by using the Box-Muller method
    {
        double a, b;
        a = rand_uniform();
        b = rand_uniopen();
        return cos(2.0*M_PI*a) * sqrt(-2.0*log(b));
    }
    
private:
    FILE *rf;
};


// ================================================================================================
bool Options::get_options(int argc, char *argv[])
{
    int c;
    bool status = true;

    while (true)
    {
        static struct option long_options[] =
        {
            {"snr", required_argument, 0, 'n'},        
            {"log2-n", required_argument, 0, 'm'},      
            {"k", required_argument, 0, 'k'},              
            {"global-multiplicity", required_argument, 0, 'M'},
            {"verbosity", required_argument, 0, 'v'},              
        };    
        
        int option_index = 0;
        c = getopt_long (argc, argv, "n:m:k:M:v:", long_options, &option_index);
        
        if (c == -1) // end of options
        {
            break;
        }

        switch(c)
        {
            case 0: // set flag
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

        rssoft::gf::GFq gfq(options.m, options.get_ppoly());
        
        URandom ur;
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
        
        // Simulate reception behind noisy channel. Reliability matrix is created with noisy power samples.

        rssoft::ReliabilityMatrix mat_Pi(options.m,n);
        std::vector<rssoft::gf::GFq_Symbol> hard_decision;
        unsigned int hard_decision_errors = 0;
        float *mat_Pi_col = new float[q];

        for (unsigned int c=0; c<n; c++)
        {
        	float max_pwr = 0;
        	unsigned int r_max = 0;

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

        delete[] mat_Pi_col;

        std::cout << "Hard-dec: (n=" << codeword.size() << ") ";
        rssoft::gf::print_symbols_vector(std::cout, hard_decision);
        std::cout << std::endl;
        std::cout << " -> " << hard_decision_errors << " errors, " << (hard_decision_errors > (n-options.k)/2 ? "uncorrectable" : "correctable") << std::endl;

        if (options.verbosity > 0)
        {
        	std::cout << "Power matrix:" << std::endl;
        	std::cout << mat_Pi;
        	std::cout << std::endl;
        }

        mat_Pi.normalize();
        float codeword_score = 0.0;
        float best_score, worst_score;

        for (unsigned int c=0; c<n; c++)
        {
        	float score = 10.0 * log10(mat_Pi(row_indexes[c],c));

        	if (c == 0)
        	{
        		best_score = score;
        		worst_score = score;
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
        }

        std::cout << "Codeword score: " << codeword_score / n << " dB/symbol (best = " << best_score << ", worst = " << worst_score << ")" << std::endl;

        rssoft::MultiplicityMatrix mat_M(mat_Pi, options.global_multiplicity);

        if (options.verbosity > 0)
        {
        	std::cout << "Multiplicity matrix:" << std::endl;
        	std::cout << mat_M;
        	std::cout << std::endl;
        }

        std::cout << "Multiplicity matrix cost is " << mat_M.cost() << std::endl;

        rssoft::GSKV_Interpolation gskv(gfq, options.k, evaluation_values);
        rssoft::RR_Factorization rr(gfq, options.k);
        gskv.set_verbosity(options.verbosity);
        rr.set_verbosity(options.verbosity);

        std::cout << std::endl;
        std::vector<rssoft::gf::GFq_Polynomial>& res_polys = rr.run(gskv.run(mat_M));

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
					std::cout << "#" << i_m << " found!!!" << std::endl;
				}
			}
        }

        return 0;
    }
    else
    {
        std::cout << "Wrong options" << std::endl;
        return -1;
    }
}
