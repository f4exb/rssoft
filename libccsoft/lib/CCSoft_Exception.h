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

     Exception for CC Soft related procedures

*/
#ifndef __CCSOFT_EXCEPTION_H__
#define __CCSOFT_EXCEPTION_H__

#include <string>

namespace ccsoft
{

    /**
     * \brief Generic exception class for Convolutional Codes Soft Decoding library.
     */
    class CCSoft_Exception : public std::exception
    {
    public:
        /**
         * Public constructor
         * \param strError Error message to be returned by the what() method
         */
        CCSoft_Exception(std::string strError) : _error_msg(strError)
        {}
        virtual ~CCSoft_Exception() throw()
        {}
        /**
         * Returns the reason for exception that was specified at exception construction time
         * \return Error message as std::string
         */
        virtual const char* what() const throw()
        {
            return _error_msg.c_str();
        }

    protected:
        std::string _error_msg; //!< Error reason shown with what() method

    private:
        /**
         *  Default constructor is not meant to be used
         */
        CCSoft_Exception()
        {};
    };
} // namespace ccsoft

#endif // __CCSOFT_EXCEPTION_H__
