lhcbSetup
SetupProject ROOT pyanalysis pytools


cflags=`root-config --cflags`
libs=`root-config --libs`

rootcflags=`root-config --cflags`
rootlibs=`root-config --libs`
allrootlibs="$ROOTLIBS -lMinuit"
# -lTreePlayer -lTMVA  -lXMLIO -lMLP -lRIO -lRooFit -lRooStats
# boost
#export LCGBOOST="$LCG_external_area/Boost/1.48.0_python2.6/$CMTCONFIG"
#export BOOST_INC_DIR="$LCGBOOST/include"
#export BOOST_LIB_OPT="-L$LCGBOOST/lib -lboost_program_options "
##-lboost_filesystem -lboost_system##-gcc46-mt-1_48 ##-lboost_system #-gcc46-mt-1_50'

# GSL
#lcggsl="$LCG_external_area/GSL/1.10/$CMTCONFIG/bin/gsl-config"
#gsldir="$lcggsl --prefix"
#gslcflags="$lcggsl --cflags"
#gslibs="$lcggsl --libs"

# libraries and flags
export LIBS=""
export CXXFLAGS=""
export CFLAGS="$rootcflags"
#$gslcflags"
export LDFLAGS="$rootlibs"
#$gsllibs"

echo $CXXFLAGS
echo $LIBS
echo $LDFLAGS
echo $CFLAGS
