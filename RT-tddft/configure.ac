#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT(reed_tddft, 1.0, sy0302@stanford.edu)

# AC_ARG_WITH( option-name, [help-string], [action-if-present], [action-if-not-present] )
# For declaring the arguments with --enable and --with prefixes, AC_ARG_ENABLE and AC_ARG_WITH work in the same way.
# AC_HELP_STRING( option-name, help-info ) gives the messaged to be displayed with ./configure --help
AC_ARG_WITH(qepath,
   [AC_HELP_STRING([--with-qepath=DIR],[Quantum Espresso source directory (QE_PATH)])],
   [qepath=$withval],
   [qepath=$QE_PATH])

# Make qe-source path absolute
# readlink -f canonicalize all symbolic links
# 2>/dev/null redirect stderr
ignore=`readlink -f ${qepath} 2>/dev/null` && qepath=$ignore

# Look for make.sys used for quantum espresso
# AC_CHECK_FILE( file, [action-if-present], [action-if-not-present] )
AC_CHECK_FILE(${qepath}/make.sys, [], [AC_MSG_ERROR([Cannot find make.sys in ${qepath}])])

# Look for quantum espresso version in qe_source/Modules/version.f90
AC_CHECK_FILE(${qepath}/Modules/version.f90, [], [AC_MSG_ERROR([Cannot find version.f90 in ${qepath}/Modules])])

# Checks for programs.

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

# Export settings to generated files
# AC_SUBST( variable, [value] )
# Make AC_OUTPUT replace instances of '@variable@' in input files with the value of the variable which can(does not have to) be reset by [value]
AC_SUBST(qepath)

# Create Makefiles from input files
# AC_CONFIG_FILES( [file] ) instruct AC_OUTPUT to create each [file] by copying an input file.in
AC_CONFIG_FILES([makedeps.sh Makefile src/Makefile])

# AC_CONFIG_COMMANDS( [tag], [cmds], [init-cmds] )
# Run [cmds] at the end of config.status with variables specified in [tag]
AC_OUTPUT
