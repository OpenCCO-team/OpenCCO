#!/bin/sh
set -e
if test "$CONFIGURATION" = "Debug"; then :
  cd /Users/kerautre/EnCours/LiverCCO/Schreiner_Algorithm/buildXcode
  make -f /Users/kerautre/EnCours/LiverCCO/Schreiner_Algorithm/buildXcode/CMakeScripts/ReRunCMake.make
fi
if test "$CONFIGURATION" = "Release"; then :
  cd /Users/kerautre/EnCours/LiverCCO/Schreiner_Algorithm/buildXcode
  make -f /Users/kerautre/EnCours/LiverCCO/Schreiner_Algorithm/buildXcode/CMakeScripts/ReRunCMake.make
fi
if test "$CONFIGURATION" = "MinSizeRel"; then :
  cd /Users/kerautre/EnCours/LiverCCO/Schreiner_Algorithm/buildXcode
  make -f /Users/kerautre/EnCours/LiverCCO/Schreiner_Algorithm/buildXcode/CMakeScripts/ReRunCMake.make
fi
if test "$CONFIGURATION" = "RelWithDebInfo"; then :
  cd /Users/kerautre/EnCours/LiverCCO/Schreiner_Algorithm/buildXcode
  make -f /Users/kerautre/EnCours/LiverCCO/Schreiner_Algorithm/buildXcode/CMakeScripts/ReRunCMake.make
fi

