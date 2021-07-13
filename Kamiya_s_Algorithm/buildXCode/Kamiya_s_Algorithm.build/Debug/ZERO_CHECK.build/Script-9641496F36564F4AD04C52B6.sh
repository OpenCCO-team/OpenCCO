#!/bin/sh
set -e
if test "$CONFIGURATION" = "Debug"; then :
  cd /Users/kerautre/EnCours/LiverCCO/Kamiya_s_Algorithm/buildXCode
  make -f /Users/kerautre/EnCours/LiverCCO/Kamiya_s_Algorithm/buildXCode/CMakeScripts/ReRunCMake.make
fi
if test "$CONFIGURATION" = "Release"; then :
  cd /Users/kerautre/EnCours/LiverCCO/Kamiya_s_Algorithm/buildXCode
  make -f /Users/kerautre/EnCours/LiverCCO/Kamiya_s_Algorithm/buildXCode/CMakeScripts/ReRunCMake.make
fi
if test "$CONFIGURATION" = "MinSizeRel"; then :
  cd /Users/kerautre/EnCours/LiverCCO/Kamiya_s_Algorithm/buildXCode
  make -f /Users/kerautre/EnCours/LiverCCO/Kamiya_s_Algorithm/buildXCode/CMakeScripts/ReRunCMake.make
fi
if test "$CONFIGURATION" = "RelWithDebInfo"; then :
  cd /Users/kerautre/EnCours/LiverCCO/Kamiya_s_Algorithm/buildXCode
  make -f /Users/kerautre/EnCours/LiverCCO/Kamiya_s_Algorithm/buildXCode/CMakeScripts/ReRunCMake.make
fi

