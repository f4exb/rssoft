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

 Debug log output macro

*/
#ifndef __DEBUG_H__
#define __DEBUG_H__

#include <time.h>

#ifdef _DEBUG
#define DEBUG_OUT(condition, str) if (condition) { std::cout << str; };
#else
#define DEBUG_OUT(condition, str) 
#endif

// time difference in seconds
double debug_get_time_difference(const timespec& time1, const timespec& time2)
{
    long long unsigned int time1_ns = time1.tv_sec * 1000000000ull + time1.tv_nsec;
    long long unsigned int time2_ns = time2.tv_sec * 1000000000ull + time2.tv_nsec;

    if (time1_ns > time2_ns)
    {
        return ((double) time1_ns - time2_ns) / 1e9;
    }
    else
    {
        return ((double) time2_ns - time1_ns) / 1e9;
    }
}

#endif // __DEBUG_H__
 
