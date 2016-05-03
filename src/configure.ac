AC_INIT([ShowerDeconstruction], [1.0], [dferreir@mail.cern.ch])
AM_INIT_AUTOMAKE([foreign subdir-objects])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES])
AC_CONFIG_MACRO_DIR([m4])

if test -z $CXXFLAGS; then
    CXXFLAGS='-O3'
fi

AC_PROG_CXX
AC_LANG(C++)

AC_PREFIX_DEFAULT(/usr/local)

AC_DEFUN([ACX_CHECK_FASTJET],
[
AC_ARG_WITH(fastjet,
            [AC_HELP_STRING([--with-fastjet=dir], 
                            [Assume the given directory for FastJet])])

if test "$with_fastjet" = ""; then
   AC_PATH_PROG(fjconfig, fastjet-config, no)
else
   AC_PATH_PROG(fjconfig, fastjet-config, no, ${with_fastjet}/bin)
fi

if test "${fjconfig}" = "no"; then
   AC_MSG_CHECKING(FastJet)
   AC_MSG_RESULT(no);
   $2
else
   dnl now see if FastJet is functional
   save_CXXFLAGS="$CXXFLAGS"
   save_LIBS="$LIBS"

   CXXFLAGS="${CXXFLAGS} `${fjconfig} --cxxflags`"
   LIBS="${LIBS} `${fjconfig} --libs --plugins`"

   AC_MSG_CHECKING([if FastJet is functional])
   AC_LANG_PUSH(C++)
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <fastjet/ClusterSequence.hh>
   ]], [[
fastjet::PseudoJet pj=fastjet::PtYPhiM(10.0,0.5,1.0,0.0);
   ]])], [fjok='yes'], [fjok='no'])
   AC_MSG_RESULT([$fjok])
   AC_LANG_POP()
   CXXFLAGS="$save_CXXFLAGS"
   LIBS="$save_LIBS"

   AC_MSG_CHECKING(FastJet)
   if test "${fjok}" = "yes"; then
      FASTJET_CXXFLAGS="`${fjconfig} --cxxflags`"
      FASTJET_LIBS="`${fjconfig} --libs --plugins`"
      AC_MSG_RESULT(yes)
      $1
   else
      AC_MSG_RESULT(no)
      $2
   fi
fi
])

ACX_CHECK_FASTJET(, AC_MSG_ERROR([Could not find a useful installation of FastJet. Try using --with-fastjet=[dir] to setup its directory properly.]))
AC_SUBST(fjlibs, ${FASTJET_LIBS})
AC_SUBST(fjcxxflags, ${FASTJET_CXXFLAGS})

AC_DISABLE_STATIC
#AC_PROG_LIBTOOL(libtool)
AC_ENABLE_SHARED
AC_PROG_LIBTOOL


AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])

AC_OUTPUT
