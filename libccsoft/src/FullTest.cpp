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

#include "ReliabilityMatrix.h"
#include "CC_Encoding.h"
#include "CC_StackDecoding.h"
#include "CC_TreeNode.h"
#include "CC_TreeEdge.h"
#include "CCSoft_Exception.h"
#include "URandom.h"

#include <getopt.h>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <sstream>
#include <cstring>

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
struct Options
{
public:
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
        generate_random_symbols(false)
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

private:
    bool parse_generator_polys_data(std::string generator_polys_data_str);
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
            // these options do not set a flag
            {"snr", required_argument, 0, 'n'},        
            {"verbosity", required_argument, 0, 'v'},              
            {"dot-output", required_argument, 0, 'd'},              
            {"k-constraints", required_argument, 0, 'k'},
            {"gen-polys", required_argument, 0, 'g'},
            {"in-symbols", required_argument, 0, 'i'},
            {"nb-random-symbols", required_argument, 0, 'r'},
            {"seed", required_argument, 0, 's'},
        };    
        
        int option_index = 0;
        c = getopt_long (argc, argv, "n:v:d:k:g:i:r:s:", long_options, &option_index);
        
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
void create_symbol_data(float *symbol_data,
        unsigned int nb_symbols,
        unsigned int out_symbol,
        float snr_dB,
        bool make_noise)
{

}

URandom ur; // Global random generator object

// ================================================================================================
int main(int argc, char *argv[])
{
    Options options;
    
    if (options.get_options(argc, argv))
    {
        try
        {
            ccsoft::CC_StackDecoding<unsigned int, unsigned int> cc_decoding(options.k_constraints, options.generator_polys);
            cc_decoding.get_encoding().print(std::cout);
            unsigned int out_symbols_nb = 1<<cc_decoding.get_encoding().get_n();
            unsigned int in_symbols_nb = 1<<cc_decoding.get_encoding().get_k();

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
                for (unsigned int i=0; i<cc_decoding.get_encoding().get_m()-1; i++)
                {
                    options.input_symbols.push_back(0);
                }
                
                ccsoft::ReliabilityMatrix relmat(cc_decoding.get_encoding().get_n(), options.input_symbols.size());
                unsigned int nb_symbols = 1<<cc_decoding.get_encoding().get_n();
                float *symbol_data = new float[nb_symbols];

                std::ostringstream oos;

                for (unsigned int i=0; i<options.input_symbols.size(); i++)
                {
                    unsigned int out_symbol;
                    cc_decoding.get_encoding().encode(options.input_symbols[i], out_symbol);
                    create_symbol_data(symbol_data, nb_symbols, out_symbol, options.snr_dB, options.make_noise);
                    relmat.enter_symbol_data(symbol_data);
                    std::cout << options.input_symbols[i] << " ";
                    oos << out_symbol << " ";
                }

                std::cout << std::endl;
                std::cout << oos.str() << std::endl;

                relmat.normalize();


                delete[] symbol_data;
            }
            
            std::cout << std::endl;
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
