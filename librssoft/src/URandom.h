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

	 Class to encapsulate /dev/urandom and rand related functions

*/


#ifndef __URANDOM_H__
#define __URANDOM_H__

#include <cmath>
#include <cstdio>

class URandom
{
public:
    URandom() : use_seed(false)
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
        
        if (use_seed)
        {
            ri = (rand()/2) - RAND_MAX;
        }
        else
        {
            unsigned int bytes_read = fread((char*)(&ri),sizeof(ri),1,rf);
        }
        
        return ri & 0x7fffffff;
    }
    
    unsigned int rand_uword() 
    {
        unsigned int ri; 
        
        if (use_seed)
        {
            ri = rand();
        }
        else
        {
            unsigned int bytes_read = fread((char*)(&ri),sizeof(ri),1,rf);
        }
        
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
    
    void set_seed(unsigned int seed)
    {
        std::cout << "use seed: " << seed << std::endl;
        srand(seed);
        use_seed = true;
    }
    
    void unset_seed()
    {
        use_seed = false;
    }
    
private:
    FILE *rf;
    bool use_seed;
};

#endif // __URANDOM_H__
