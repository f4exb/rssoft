AC_INIT(rssoft, 0.0.0)
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE(rssoft, 0.0.0)
LT_INIT
AC_LANG_PUSH([C++])

ARCH=`uname -m`

# ------------------------------------------------------------------------------
# DOXYGEN
# ------------------------------------------------------------------------------

AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN";
   then AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])
fi

AM_CONDITIONAL([HAVE_DOXYGEN], 
[test -n "$DOXYGEN"])AM_COND_IF([HAVE_DOXYGEN], [AC_CONFIG_FILES([doc/Doxyfile])])

dnl: this is a comment
dnl: ${CXXFLAGS=-g}
AC_PROG_CXX

# ------------------------------------------------------------------------------
# NTL library
# ------------------------------------------------------------------------------
    
dnl:AC_ARG_WITH(ntl,
dnl:  [AS_HELP_STRING([--with-ntl=<path>],
dnl:                  [prefix of NTL installation. e.g. /usr/local or /usr])],
dnl:  [],
dnl:  [AC_CHECK_LIB([ntl], [_ntl_gmul],[ ],
dnl:    AC_MSG_ERROR([NTL Library not found.
dnl:    Maybe you want to call configure with the --with-ntl=<path> option?
dnl:    This tells configure where to find the NTL library and headers.
dnl:    e.g. --with-ntl=/usr/local or --with-ntl=/usr]),[${GMP_LIBS}])]
dnl:)
dnl:case "$with_ntl" in
dnl:  ""|yes)
dnl:    NTL_LDFLAGS=""
dnl:    NTL_CFLAGS=""
dnl:    ;;
dnl:  no)
dnl:    AC_MSG_ERROR([Compiling without libntl impossible (but --without-ntl given)])
dnl:    ;;
dnl:  *)
dnl:    NTL_LDFLAGS="-L${with_ntl}/lib"
dnl:    NTL_CFLAGS="-I${with_ntl}/include"
dnl:    ;;
dnl:esac
dnl:NTL_LIBS="-lntl"
dnl:AC_SUBST(NTL_LIBS)
dnl:AC_SUBST(NTL_CFLAGS)
dnl:AC_SUBST(NTL_LDFLAGS)

# ------------------------------------------------------------------------------
# EPILOG
# ------------------------------------------------------------------------------
    
AC_CONFIG_FILES([Makefile doc/Makefile lib/Makefile src/Makefile])
AC_OUTPUT


