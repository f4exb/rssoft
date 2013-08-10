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
#include <getopt.h>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

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
// template to extract a vector of elements from a comma separated string
template<typename TElement> bool extract_vector(std::vector<TElement>& velements, std::string cs_string)
{
    std::string element_str;
    TElement element;

    boost::char_separator<char> sep(",");
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
        std::cout << "wrong element in comma separated string argument: " << *tok_iter << std::endl;
        return false;
    }
}

// ================================================================================================
int main(int argc, char *argv[])
{
    return 0;
}
