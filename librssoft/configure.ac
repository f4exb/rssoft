AC_INIT(rssoft, 0.0.0)
AC_CONFIG_MACRO_DIR([m4])
m4_include([m4/boost.m4])
AM_INIT_AUTOMAKE(rssoft, 0.0.0)
LT_INIT
AC_LANG_PUSH([C++])

ARCH=`uname -m`

AC_PROG_CXX

# ------------------------------------------------------------------------------
# BOOST
# ------------------------------------------------------------------------------
BOOST_REQUIRE([1.33])

# ------------------------------------------------------------------------------
# EPILOG
# ------------------------------------------------------------------------------
    
dnl: this is a comment
dnl: ${CXXFLAGS=-g}

AC_CONFIG_FILES([Makefile lib/Makefile src/Makefile])
AC_OUTPUT


