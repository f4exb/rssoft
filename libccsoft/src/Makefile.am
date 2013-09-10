AM_CPPFLAGS = -I$(srcdir)/../lib
bin_PROGRAMS = Encoder_test Decoder_test FullTest FullTest_FA Sizes Interleaver_test

Encoder_test_SOURCES = Encoder_test.cpp
Encoder_test_LDADD = ../lib/libccsoft.la

Decoder_test_SOURCES = Decoder_test.cpp
Decoder_test_LDADD = ../lib/libccsoft.la

Interleaver_test_SOURCES = Interleaver_test.cpp
Interleaver_test_LDADD = ../lib/libccsoft.la

FullTest_SOURCES = FullTest.cpp
FullTest_CPPFLAGS = -I$(srcdir)/../lib $(BOOST_CPPFLAGS)
FullTest_LDADD = ../lib/libccsoft.la -lrt

FullTest_FA_SOURCES = FullTest_FA.cpp
FullTest_FA_CPPFLAGS = -std=c++0x -I$(srcdir)/../lib $(BOOST_CPPFLAGS)
FullTest_FA_LDADD = ../lib/libccsoft.la -lrt

Sizes_SOURCES = Sizes.cpp
Sizes_CPPFLAGS = -std=c++0x -I$(srcdir)/../lib $(BOOST_CPPFLAGS)
Sizes_LDADD = ../lib/libccsoft.la -lrt
