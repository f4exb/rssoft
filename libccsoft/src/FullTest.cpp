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

#include <getopt.h>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <iostream>

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
        nb_random_symbols(0)
    {}
    
    ~Options()
    {}
    
    bool get_options(int argc, char *argv[]);
        
    bool make_noise;
    bool dot_output;
    float snr_dB;
    unsigned int verbosity;
    unsigned int nb_random_symbols;
    std::string dot_filename;
    std::vector<unsigned int> k_constraints;
    std::vector<std::vector<unsigned int> > generator_polys;
    std::vector<unsigned int> input_symbols;
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
            // these options do not set a flag
            {"snr", required_argument, 0, 'n'},        
            {"verbosity", required_argument, 0, 'v'},              
            {"dot-output", required_argument, 0, 'd'},              
            {"k-constraints", required_argument, 0, 'k'},
            {"gen-polys", required_argument, 0, 'g'},
            {"in-symbols", required_argument, 0, 'i'},
            {"nb-random-symbols", required_argument, 0, 'r'},
        };    
        
        int option_index = 0;
        c = getopt_long (argc, argv, "n:v:d:k:g:i:r:", long_options, &option_index);
        
        if (c == -1) // end of options
        {
            break;
        }

        switch(c)
        {
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
int main(int argc, char *argv[])
{
    Options options;
    
    if (options.get_options(argc, argv))
    {
        try
        {
            ccsoft::CC_StackDecoding<unsigned int, unsigned int> cc_decoding(options.k_constraints, options.generator_polys);
            cc_decoding.get_encoding().print(std::cout);
            
            if (options.input_symbols.size() > 0)
            {
                for (unsigned int i=0; i<cc_decoding.get_encoding().get_m()-1; i++)
                {
                    options.input_symbols.push_back(0);
                }
                
                for (unsigned int i=0; i<options.input_symbols.size(); i++)
                {
                    unsigned int out_symbol;
                    cc_decoding.get_encoding().encode(options.input_symbols[i], out_symbol);
                    std::cout << out_symbol << " ";
                }
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
