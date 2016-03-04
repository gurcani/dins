#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)

AC_INIT([dins], [0.1], [ozgur.gurcan@lpp.polytechnique.fr])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_CONFIG_SRCDIR([])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_RANLIB

AC_ARG_ENABLE(debug, AC_HELP_STRING([--enable-debug],[enables debugging]), ok=$enableval, ok=no)
if test "$ok" = "yes"; then
        AC_DEFINE(SHELLPP_DEBUG,1,[define to enable debugging code.])
        CXXFLAGS=" -g -O0"
else
        CXXFLAGS=" -O3"
fi


# Checks for libraries.
AC_CHECK_LIB([hdf5],[main],[CXXFLAGS+=" -DH5_USE_16_API" LIBS+="-lhdf5 " ],[echo "Error! You need hdf5."
                        echo "see http://hdf.ncsa.uiuc.edu/HDF5/"
                        exit -1
                        ])

AC_CHECK_LIB([hdf5_hl],[main],[],[echo "Error! You need hdf5_hl."
                        echo "see http://hdf.ncsa.uiuc.edu/HDF5/"
                        exit -1
                        ])

AC_CHECK_LIB([boost_system],[main],[],[echo "Error! You need boost."
                        echo "see http://www.boost.org/"
                        exit -1
                        ])
AC_CHECK_LIB([boost_timer],[main],[],[echo "Error! You need boost timer libs."
                        echo "see http://www.boost.org/doc/libs/1_59_0/libs/timer/doc/index.html"
                        exit -1
                        ])

#AC_CHECK_LIB([m],[main],[],[echo "Error! You need libm."
#                        echo  "see www.gnu.org/software/libc/"
#                        exit -1
#                        ])

#AC_CHECK_LIB([gslcblas],[main],[LIBS+="-lgsl -lgslcblas"],[echo "Error! You need gsl."
#                        echo  "see http://www.gnu.org/software/gsl/"
#                        exit -1
#                        ])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([boost/array.hpp])
AC_CHECK_HEADERS([boost/numeric/odeint.hpp])
AC_CHECK_HEADERS([omp.h],[CXXFLAGS+=" -fopenmp"],[])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
#AC_FUNC_MALLOC

#AC_CONFIG_FILES([Makefile])
AC_OUTPUT(Makefile)
