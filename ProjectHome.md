

# Introduction #
List decoding and soft decision based algorithms for the Reed-Solomon codes were introduced by the works of Madhu Sudan, Venkatesan Guruswami, Ralf Koetter and Alexander Vardy.

There were a few academic papers introducing their work in more details and more noticeably the [thesis](http://sist.sysu.edu.cn/~chenli/Thesis%20(Li%20Chen).pdf) of Li Chen who proposed an optimisation to the Guruswami-Sudan interpolation algorithm.

# Documentation #
## Reed-Solomon codes part 1: introduction, classical and generic decoding ##
This part in [English](https://code.google.com/p/rssoft/downloads/detail?name=rs_part1_classical_decoding.en.pdf) and [French](https://code.google.com/p/rssoft/downloads/detail?name=rs_part1_classical_decoding.fr.pdf) versions introduces the basic mathematical notions in use with Reed-Solomon codes. It then explains the coding and a classical method of decoding based on syndrome calculation. Eventually a less classical generic decoding is presented that can be used with both systematic and non-systematic Reed-Solomon codes.
## Reed-Solomon codes part 2: soft decision decoding ##
This part in [English](https://code.google.com/p/rssoft/downloads/detail?name=rs_part2_soft_decoding.en.pdf) and [French](https://code.google.com/p/rssoft/downloads/detail?name=rs_part2_soft_decoding.fr.pdf) versions aims at explaining how the soft decision based decoding of Reed-Solomon codes works based on various sources and compiled into a synthetic, simplified and hopefully useful document. This is the base of the software implementation hosted here.
## Library Doxygen Documentation ##
Library documentation is accessible [here](http://f4exb.free.fr/rssoft/library/doc/html/index.html)

# Development status #
I haven't started the "official" releasing of the library and therefore it is at version 0.0.0 according to libtool scheme which means draft status. Although it has the full functionnality I am waiting to exercise it in a project closer to the real life of weak signal processing (see "wsgc" project in Google Code).

At the moment it is exclusively targeted to Linux systems. Feel free to contact me if you want to collaborate for making it available to other environments such as MS-Windows. This is plain C++ code and the library code itself does not adhere to any other library except the standard library (STL) so that shouldn't be too much of a problem.

**New!** A library for convolutional codes soft-decision decoding is being implemented as complementary to the original Reed-Solomon one. Often both types are necessary and it could be nice to have both in one place but in different libraries. Check the "libccsoft" directory in the source tree.

# Source contents #
## Sage scripts ##
A few Sage scripts demonstrate the essential steps of soft-decoding Reed-Solomon codes and serve as reference for unit testing the C++ library:
  * Guruswami-Sudan modified by Koetter-Vardy interpolation algorithm with Li Chen's optimization
  * Roth-Ruckenstein factorization algorithm with my own optimization in the pure computational domain

## Libraries ##
### Reed-Solomon library ###
The "librssoft' sub-directory comprises:
  * A complete Galois Field layer for our purpose. As I am not quite satisfied with the code style and/or the functionality of existing libraries, I have written my own inspired by several sources. It includes support for bivariate polynomials as required by the various algorithms.
  * Reliability matrix and Multiplicity matrix construction
  * The Guruswami-Sudan-Koetter-Vardy interpolation
  * The Roth-Ruckenstein factorization
  * The final evaluation of possible codewords and messages
All pass unit testing (gives same results as Sage prototypes).
To see a complete example of usage you may look at the [decoding unit test](https://code.google.com/p/rssoft/source/browse/library/src/Decode_UnitTest.cpp)

Doxygen documentation has been posted [here](http://f4exb.free.fr/rssoft/library/doc/html/index.html)

### Convolutional codes library ###
The "libccsoft" sub-directory contains an encoder and various flavours of decoders. They are templated classes capable of supporting various integer types for internal registers and input/output symbols. It supports multiple input lines (k>1) or input symbol size and output lines (n) or output symbol size of course as long as k<n. The constraint length is also parametrizable through the internal register lengths.

Decoding has two flavours able to handle long constraints. Therefore there is no Viterbi implementation.
  * Stack decoding: based on the stack algorithm. It is mainly memory constrained as the stack and code tree size grows exponentially large as we approach the limit SNR for successful decoding.
  * Fano decoding: based on the Fano algorithm. It is mainly processor constrained. There is a minimal footprint as only the current path is kept in memory at the expense of more processing with constant edge and nodes construction and destruction. Nodes and edges can be optionnally kept in memory by specifying a maximum number of cached nodes. When the limit is reached and new nodes and edges need to be allocated memory is purged of all nodes and edges except those on the current path and the immediate successors of all nodes in the path.

Doxygen documentation has been posted [here](http://f4exb.free.fr/ccsoft/library/doc/html/index.html)

## Test programs ##
### Reed-Solomon ###
  * Full test program: "FullTest" executable gets installed in the bin sub-directory. For details see [wiki page here](https://code.google.com/p/rssoft/wiki/FullTest_program)
  * Unit test using the Sage scripts example. see [here](https://code.google.com/p/rssoft/source/browse/library/src/Decode_UnitTest.cpp)

### Convolutional ###
  * Encoder unit test: Encoder\_test
  * Decoder unit test: Decoder\_test
  * Full test: [documentation here](https://code.google.com/p/rssoft/wiki/FullTest_ccsoft)
  * Convolutional codes tables:
    * [http://ipnpr.jpl.nasa.gov/progress\_report/42-80/80K.PDF](http://ipnpr.jpl.nasa.gov/progress_report/42-80/80K.PDF)

## Running the Reed-Solomon full test program ##

The solution appears to run globally as expected however there are a few disturbing facts. Read [this article in the wiki](https://code.google.com/p/rssoft/wiki/AFewObservations)

## Limitations of the Reed-Solomon decoder ##

The current code can only work with generalized Reed-Solomon codes (GRS) not systematic codes.

# Comparison with Koetter and Vardy's KV (a.k.a. KVASD) algorithm #

See the test done with WSGC simulator in WSGC wiki [here](https://code.google.com/p/wsgc/wiki/RSSoft_vs_KV)
Results show that performance is worse by about 0.5 dB.