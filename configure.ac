#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([ga], [1.3], [yosh36@outlook.com])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE


# Checks for programs.
AC_PROG_MAKE_SET
AC_PROG_CXX
AC_LANG(C++)

CXXFLAGS="-O2 -std=c++11 -lboost_timer"

AC_ENABLE_SHARED
AC_CONFIG_MACRO_DIR([m4])
LT_PREREQ([2.2])
LT_INIT([dlopen shared])

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt])

AC_CONFIG_FILES([Makefile src/Makefile test/Makefile])
AC_OUTPUT
