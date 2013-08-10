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

 Reed-Solomon encoding class

 */
#include "RS_SystematicEncoding.h"
#include "GFq.h"
#include "RSSoft_Exception.h"

namespace rssoft
{

// ================================================================================================
RS_SystematicEncoding::RS_SystematicEncoding(const gf::GFq& _gf, unsigned int _k, unsigned int _init_power) :
	gf(_gf),
	k(_k),
	init_power(_init_power),
	G(_gf)
{
	// Construct generator polynomial

	// X-a^i
	std::vector<gf::GFq_Element> xe;
	xe.push_back(gf::GFq_Element(gf,gf.alpha(init_power)));
	xe.push_back(gf::GFq_Element(gf,1));

	G.init(xe);
	gf::GFq_Polynomial X(gf);

	for (unsigned int i=1; i < gf.size() - k; i++)
	{
		xe[0] = gf::GFq_Element(gf,gf.alpha(init_power+i));
		X.init(xe);
		G *= X;
	}
}

// ================================================================================================
RS_SystematicEncoding::~RS_SystematicEncoding()
{}

// ================================================================================================
void RS_SystematicEncoding::run(const std::vector<gf::GFq_Symbol>& message, std::vector<gf::GFq_Symbol>& codeword) const
{
	if (message.size() != k)
	{
		throw RSSoft_Exception("Invalid message length");
	}
	else
	{
		std::cout << "G(X) = " << G << std::endl;
		std::vector<gf::GFq_Symbol>::const_iterator s_it = message.begin();
		std::vector<gf::GFq_Element> encoding_coefficients;

		for(; s_it != message.end(); ++s_it)
		{
			encoding_coefficients.push_back(gf::GFq_Element(gf, *s_it));
		}

		gf::GFq_Polynomial encoding_polynomial(gf, encoding_coefficients);
		gf::GFq_Polynomial X_nk(gf::GFq_Element(gf,1), gf.size() - k);
		gf::GFq_Polynomial shifted_encoding_polynomial = X_nk*encoding_polynomial;
		std::pair<gf::GFq_Polynomial, gf::GFq_Polynomial> Q_R = gf::div(shifted_encoding_polynomial, G);
		gf::GFq_Polynomial codeword_polynomial = Q_R.second + shifted_encoding_polynomial;
		codeword_polynomial.get_poly_symbols(codeword);
	}
}


// ================================================================================================
/*
void RS_SystematicEncoding::encode_rs()
// Compute the 2t parity symbols in b[0]..b[2*t-1]
// data[] is input and b[] is output in polynomial form.
// Encoding is done by using a feedback shift register with connections
// specified by the elements of g[].
{
   register int i,j;
   unsigned int feedback;
   unsigned int length = gf.size();
   unsigned int recd[63], data[12], b[51];

   for (i=0; i<length-k; i++)
     b[i] = 0;
   for (i=k-1; i>=0; i--)
    {
    //feedback = index_of[data[i]^b[length-k-1]];
    feedback = gf.index(data[i]^b[length-k-1]);
      if (feedback != -1)
        {
        for (j=length-k-1; j>0; j--)
          if (g[j] != -1)
            b[j] = b[j-1]^gf.alpha((g[j]+feedback)%n);
          else
            b[j] = b[j-1];
          b[0] = alpha_to[(g[0]+feedback)%n];
        }
       else
        {
        for (j=length-k-1; j>0; j--)
          b[j] = b[j-1];
        b[0] = 0;
        }
    }
}
*/

} // namespace rssoft
