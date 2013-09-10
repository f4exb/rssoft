/*
     Copyright 2013 Edouard Griffiths <f4exb at free dot fr>

     This file is part of CCSoft. A Convolutional Codes Soft Decoding library

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

#include "CC_ReliabilityMatrix.h"
#include "CC_Encoding.h"
#include "CC_StackDecoding.h"
#include "CC_FanoDecoding.h"
#include "CCSoft_Exception.h"
#include "URandom.h"

#include <getopt.h>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <cctype> // for toupper

static URandom ur; // Global random generator object


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
// template to extract a vector of elements from a delimiter separated string
template<typename TElement> bool extract_vector(std::vector<TElement>& velements, const char *separator, std::string cs_string)
{
    std::string element_str;
    TElement element;

    boost::char_separator<char> sep(separator);
    boost::tokenizer<boost::char_separator<char> > tokens(cs_string, sep);

    boost::tokenizer<boost::char_separator<char> >::iterator tok_iter = tokens.begin();
    boost::tokenizer<boost::char_separator<char> >::iterator toks_end = tokens.end();

    try
    {
        for (; tok_iter != toks_end; ++tok_iter)
        {
            element = boost::lexical_cast<TElement>(*tok_iter);
            velements.push_back(element);
        }
        return true;
    }
    catch (boost::bad_lexical_cast &)
    {
        std::cout << "wrong element in delimiter separated string argument: " << *tok_iter << std::endl;
        return false;
    }
}

// ================================================================================================
// template to print a vector of printable elements
template<typename TDisplay, typename TElement, typename TStream> void print_vector(const std::vector<TElement>& v, TStream& os)
{
    os << "[";

    typename std::vector<TElement>::const_iterator it = v.begin();
    const typename std::vector<TElement>::const_iterator v_end = v.end();

    for (; it != v_end; ++it)
    {

        os <<  (TDisplay)(*it);

        if (it != v.begin()+v.size()-1)
        {
            os << ", ";
        }

    }

    os << "]";
}


// ================================================================================================
struct Options
{
public:
	typedef enum
	{
		Algorithm_Stack,
		Algorithm_FanoLike
	} Algorithm_type_t;

    Options() :
        make_noise(false),
        dot_output(false),
        snr_dB(0),
        verbosity(0),
        indicator_int(0),
        print_seed(false),
        seed(0),
        has_seed(false),
        nb_random_symbols(0),
        generate_random_symbols(false),
        node_limit(0),
        use_node_limit(false),
        metric_limit(0.0),
        use_metric_limit(false),
        algorithm_type(Algorithm_Stack),
        fano_init_metric(-1.0),
        fano_delta_metric(1.0),
        fano_tree_cache_size(0),
        edge_bias(0.0),
        fano_delta_init_threshold(0.0),
        interleave(false)
    {}

    ~Options()
    {}

    bool get_options(int argc, char *argv[]);

    bool make_noise;
    bool dot_output;
    float snr_dB;
    unsigned int verbosity;
    std::string dot_filename;
    std::vector<unsigned int> k_constraints;
    std::vector<std::vector<unsigned int> > generator_polys;
    std::vector<unsigned int> input_symbols;
    int indicator_int;
    bool print_seed;
    unsigned int seed;
    bool has_seed;
    unsigned int nb_random_symbols;
    bool generate_random_symbols;
    unsigned int node_limit;
    bool use_node_limit;
    float metric_limit;
    bool use_metric_limit;
    Algorithm_type_t algorithm_type;
    float fano_init_metric;
    float fano_delta_metric;
    unsigned int fano_tree_cache_size;
    float edge_bias;
    float fano_delta_init_threshold;
    bool interleave;

private:
    bool parse_generator_polys_data(std::string generator_polys_data_str);
    bool parse_algorithm_type(std::string algorithm_type_str);
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
            // these options set a flag
            {"print-seed", no_argument, &indicator_int, 1},
            {"interleave", no_argument, &indicator_int, 1},
            // these options do not set a flag
            {"snr", required_argument, 0, 'n'},
            {"verbosity", required_argument, 0, 'v'},
            {"dot-output", required_argument, 0, 'd'},
            {"k-constraints", required_argument, 0, 'k'},
            {"gen-polys", required_argument, 0, 'g'},
            {"in-symbols", required_argument, 0, 'i'},
            {"nb-random-symbols", required_argument, 0, 'r'},
            {"seed", required_argument, 0, 's'},
            {"node-limit", required_argument, 0, 'N'},
            {"metric-limit", required_argument, 0, 'M'},
            {"algorithm-type", required_argument,0, 'a'},
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "n:v:d:k:g:i:r:s:N:M:a:", long_options, &option_index);

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
                else if (strcmp("interleave", long_options[option_index].name) == 0)
                {
                    interleave = true;
                }
                break;
            case 'n':
                make_noise = true;
                status = extract_option<double, float>(snr_dB, 'n');
                break;
            case 'v':
                status = extract_option<int, unsigned int>(verbosity, 'v');
                break;
            case 'd':
                dot_output = true;
                dot_filename = std::string(optarg);
                break;
            case 'k':
                status = extract_vector<unsigned int>(k_constraints, ",", std::string(optarg));
                break;
            case 'g':
                status = parse_generator_polys_data(std::string(optarg));
                break;
            case 'i':
                status = extract_vector<unsigned int>(input_symbols, ",", std::string(optarg));
                break;
            case 'r':
                status = extract_option<int, unsigned int>(nb_random_symbols, 'r');
                generate_random_symbols = true;
                break;
            case 's':
                status = extract_option<int, unsigned int>(seed, 's');
                has_seed = true;
                break;
            case 'N':
                status = extract_option<int, unsigned int>(node_limit, 'N');
                use_node_limit = true;
                break;
            case 'M':
                status = extract_option<float, float>(metric_limit, 'M');
                use_metric_limit = true;
                break;
            case 'a':
                status = parse_algorithm_type(std::string(optarg));
                break;
            case '?':
                status = false;
                break;
        }
    }
}

// ================================================================================================
bool Options::parse_generator_polys_data(std::string generator_polys_data_str)
{
    std::vector<std::string> g_strings;

    if (!extract_vector(g_strings, ":", generator_polys_data_str))
    {
    	std::cerr << "Invalid generator polynomials specification" << std::endl;
        return false;
    }

    std::vector<std::string>::const_iterator gs_it = g_strings.begin();

    for (; gs_it != g_strings.end(); ++gs_it)
    {
        std::vector<unsigned int> g;

        if (extract_vector<unsigned int>(g, ",", *gs_it))
        {
            generator_polys.push_back(g);
        }
        else
        {
            return false;
        }
    }

    return true;
}

// ================================================================================================
bool Options::parse_algorithm_type(std::string algorithm_type_str)
{
    std::vector<std::string> algo_strings;

    if (!extract_vector(algo_strings, ":", algorithm_type_str))
    {
    	std::cerr << "Invalid algorithm specification" << std::endl;
        return false;
    }

    std::transform(algo_strings[0].begin(), algo_strings[0].end(), algo_strings[0].begin(), toupper);

	if (algo_strings[0] == "FANO")
	{
		if (algo_strings.size() > 1)
		{
			std::vector<float> fano_parms;

			if (extract_vector(fano_parms, ",", algo_strings[1]))
			{
				if (fano_parms.size() > 0)
				{
					edge_bias = fano_parms[0];
				}
				if (fano_parms.size() > 1)
				{
					fano_init_metric = fano_parms[1];
				}
				if (fano_parms.size() > 2)
				{
					fano_delta_metric = fano_parms[2];
				}
                if (fano_parms.size() > 3)
                {
                    fano_tree_cache_size = int(fano_parms[3]);
                }
                if (fano_parms.size() > 4)
                {
                    fano_delta_init_threshold = int(fano_parms[4]);
                }
			}
			else
			{
				std::cerr << "Invalid Fano parameters specification" << std::endl;
				return false;
			}
		}

		algorithm_type = Algorithm_FanoLike;
		return true;
	}
	else if (algo_strings[0] == "STACK")
	{
		if (algo_strings.size() > 1)
		{
			std::vector<float> stack_parms;

			if (extract_vector(stack_parms, ",", algo_strings[1]))
			{
				if (stack_parms.size() > 0)
				{
					edge_bias = stack_parms[0];
				}
			}
			else
			{
				std::cerr << "Invalid Stack parameters specification" << std::endl;
				return false;
			}
		}

		algorithm_type = Algorithm_Stack;
		return true;
	}
	else
	{
		return false;
	}
}

// ================================================================================================
void create_symbol_data(float *symbol_data,
        unsigned int nb_symbols,
        unsigned int out_symbol,
        float snr_dB,
        bool make_noise)
{
    double std_dev  = 1.0 / pow(10.0, (snr_dB/10.0)); // Standard deviation for power AWGN

    for (unsigned int si=0; si<nb_symbols; si++)
    {
        if (si == out_symbol)
        {
            symbol_data[si] = 1.0 + (make_noise ? std_dev * ur.rand_gaussian() : 0.0);
        }
        else
        {
            symbol_data[si] = 0.0 + (make_noise ? std_dev * ur.rand_gaussian() : 0.0);
        }

        symbol_data[si] *= symbol_data[si];
    }
}

// ================================================================================================
int main(int argc, char *argv[])
{
    Options options;
    bool success = false;

    if (options.get_options(argc, argv))
    {
        ccsoft::CC_SequentialDecoding<unsigned int, unsigned int> *cc_decoding;
    
        try
        {
            if (options.algorithm_type == Options::Algorithm_Stack)
            {
                cc_decoding = new ccsoft::CC_StackDecoding<unsigned int, unsigned int>(options.k_constraints, options.generator_polys);
            }
            else if (options.algorithm_type == Options::Algorithm_FanoLike)
            {
                cc_decoding = new ccsoft::CC_FanoDecoding<unsigned int, unsigned int>(options.k_constraints,
                        options.generator_polys,
                        options.fano_init_metric,
                        options.fano_delta_metric,
                        options.fano_tree_cache_size,
                        options.fano_delta_init_threshold);
            }
            else
            {
                std::cerr << "Unrecognized algorithm type" << std::endl;
                return 1;
            }
             
            cc_decoding->set_verbosity(options.verbosity);
            cc_decoding->set_edge_bias(options.edge_bias); 
            cc_decoding->get_encoding().print(std::cout);
            unsigned int out_symbols_nb = 1<<cc_decoding->get_encoding().get_n();
            unsigned int in_symbols_nb = 1<<cc_decoding->get_encoding().get_k();

            if (options.use_node_limit)
            {
                cc_decoding->set_node_limit(options.node_limit);
            }

            if (options.use_metric_limit)
            {
                cc_decoding->set_metric_limit(options.metric_limit);
            }

            if (options.has_seed)
            {
                ur.set_seed(options.seed);
            }

            if (options.print_seed)
            {
                unsigned int seed = ur.rand_uword();
                std::cout << "Seed = " << std::dec << seed << std::endl;
                ur.set_seed(seed);
            }

            if (options.generate_random_symbols)
            {
                options.input_symbols.clear(); // ignore given input symbols

                for (unsigned int i=0; i<options.nb_random_symbols; i++)
                {
                    options.input_symbols.push_back(ur.rand_int(in_symbols_nb));
                }
            }

            if (options.input_symbols.size() > 0)
            {
                for (unsigned int i=0; i<cc_decoding->get_encoding().get_m()-1; i++)
                {
                    options.input_symbols.push_back(0);
                }

                ccsoft::CC_ReliabilityMatrix relmat(cc_decoding->get_encoding().get_n(), options.input_symbols.size());
                unsigned int nb_symbols = 1<<cc_decoding->get_encoding().get_n();
                float *symbol_data = new float[nb_symbols];

                std::ostringstream oos;

                if (options.interleave)
                {
                	std::cout << "interleave" << std::endl;
                	std::vector<unsigned int> out_symbols;

					for (unsigned int i=0; i<options.input_symbols.size(); i++)
					{
						unsigned int out_symbol;
						cc_decoding->get_encoding().encode(options.input_symbols[i], out_symbol);
						out_symbols.push_back(out_symbol);
						std::cout << options.input_symbols[i] << " ";
						oos << out_symbol << " ";
					}

					cc_decoding->interleave(out_symbols);

					for (unsigned int i=0; i<out_symbols.size(); i++)
					{
						create_symbol_data(symbol_data, nb_symbols, out_symbols[i], options.snr_dB, options.make_noise);
						relmat.enter_symbol_data(symbol_data);
					}

					relmat.deinterleave();
                }
                else
                {
					for (unsigned int i=0; i<options.input_symbols.size(); i++)
					{
						unsigned int out_symbol;
						cc_decoding->get_encoding().encode(options.input_symbols[i], out_symbol);
						create_symbol_data(symbol_data, nb_symbols, out_symbol, options.snr_dB, options.make_noise);
						relmat.enter_symbol_data(symbol_data);
						std::cout << options.input_symbols[i] << " ";
						oos << out_symbol << " ";
					}
                }

                delete[] symbol_data;
                std::cout << std::endl;
                std::cout << oos.str() << std::endl;

                relmat.normalize();
                std::vector<unsigned int> result;

                if (cc_decoding->decode(relmat, result))
                {
                    print_vector<unsigned int>(result, std::cout);
                    std::cout << " ";

                    success = (result == options.input_symbols);
                    
                    if (success)
                    {
                        std::cout << "Success!" << std::endl;
                    }
                    else
                    {
                        std::cout << "Failed :(" << std::endl;
                    }

                    if (options.dot_output)
                    {
                        std::ofstream dot_file;
                        dot_file.open(options.dot_filename.c_str());
                        cc_decoding->print_dot(dot_file);
                        dot_file.close();
                    }
                }
                else
                {
                    std::cout << "Message cannot be decoded" << std::endl;
                }
                
                cc_decoding->print_stats(std::cout, success);
            }
        }
        catch (ccsoft::CCSoft_Exception& e)
        {
            std::cout << "CCSoft exception caught: " << e.what() << std::endl;
        }
    }
    else
    {
        std::cout << "Wrong options" << std::endl;
        return -1;
    }

    return 0;
}
