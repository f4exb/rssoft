AM_CPPFLAGS = -I$(srcdir)/../lib
bin_PROGRAMS = GF8_test GF2_test GF8_bpoly_test Decode_UnitTest FullTest

GF8_test_SOURCES = GF8_test.cpp
GF8_test_LDADD = ../lib/librssoft.la

GF2_test_SOURCES = GF2_test.cpp
GF2_test_LDADD = ../lib/librssoft.la

GF8_bpoly_test_SOURCES = GF8_bpoly_test.cpp
GF8_bpoly_test_LDADD = ../lib/librssoft.la

Decode_UnitTest_SOURCES = Decode_UnitTest.cpp
Decode_UnitTest_LDADD = ../lib/librssoft.la

FullTest_SOURCES = FullTest.cpp
FullTest_CPPFLAGS = -I$(srcdir)/../lib $(BOOST_CPPFLAGS)
FullTest_LDADD = ../lib/librssoft.la
