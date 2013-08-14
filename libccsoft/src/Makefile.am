AM_CPPFLAGS = -I$(srcdir)/../lib
bin_PROGRAMS = Encoder_test Decoder_test FullTest

Encoder_test_SOURCES = Encoder_test.cpp
Encoder_test_LDADD = ../lib/libccsoft.la

Decoder_test_SOURCES = Decoder_test.cpp
Decoder_test_LDADD = ../lib/libccsoft.la

FullTest_SOURCES = FullTest.cpp
FullTest_CPPFLAGS = -I$(srcdir)/../lib $(BOOST_CPPFLAGS)
FullTest_LDADD = ../lib/libccsoft.la
