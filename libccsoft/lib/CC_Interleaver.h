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

 Interleaver class

 */
#ifndef __CC_INTERLEAVER_H__
#define __CC_INTERLEAVER_H__

#include <cmath>
#include <iostream>

namespace ccsoft
{

template<typename T_IOSymbol>
class CC_Interleaver
{
public:
    /**
     * Interleave/Deinterleave
     * \param symbols Symbols to process
     * \param true if forward direction i.e. interleave, false for de-interleave
     */
    void interleave(std::vector<T_IOSymbol>& symbols, bool forward=true)
    {
        std::vector<T_IOSymbol> tmp_symbols(symbols);
        unsigned int index_size = (unsigned int) (log(symbols.size())/log(2)) + 1;
        unsigned int index_max = 1<<index_size;
        unsigned int new_index, s, iv;
        unsigned int old_index = 0;

        for (unsigned int i=0; (i<index_max) && (old_index<symbols.size()); i++)
        {
            new_index = 0;
            s = index_size;
            iv = i;

            for (; iv; iv >>= 1) // bit reversal
            {
                new_index |= iv & 1;
                new_index <<= 1;
                s--;
            }

            new_index >>= 1; // the last shift right was too much
            new_index <<= s; // account for leading zeroes

            if (new_index < symbols.size())
            {
                if (forward)
                {
                    symbols[new_index] = tmp_symbols[old_index];
                }
                else
                {
                	symbols[old_index] = tmp_symbols[new_index];
                }

                old_index++;
            }
        }
    }
};

} // namespace ccsoft

#endif // __CC_INTERLEAVER_H__
