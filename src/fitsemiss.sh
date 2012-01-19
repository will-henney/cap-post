##
## Make emissivity datacube FITS files
##
## 01 Dec 2005 WJH
##
## Example usage:
##
##        sh fitsemiss.sh 30112005 c Halpha
##

dateid=$1
runid=$2
emtype=$3

# optionally restrict lowest and highest time
if [ $# -ge 4 ]; then
    i1=$4
else
    i1=1
fi
if [ $# -ge 5 ]; then
    i2=$5
else
    i2=9999
fi
if [ $# -ge 6 ]; then
    iskip=$6
else
    iskip=1
fi


## This is simpler than multifits.sh since we assume all data
## files are in the current directory
prefix1=$(printf '%s_%s' $dateid $runid)

execdir=$(dirname $0)

for i in $(seq $i1 $iskip $i2); do
    prefix2=$(printf '%s_%4.4i' $prefix1 $i)
    if [ -f ${prefix2}"d.fits" ]; then
	if [ ! -f ${prefix2}"t.fits" ]; then
	    echo "Calculating temperature cube for $prefix2"
	    printf '%s\n' $prefix2 | ${execdir}/cubet > /dev/null 
	fi
	if [ ! -f ${prefix2}"e-"${emtype}".fits" ]; then
	    echo "Calculating $emtype emissivity cube for $prefix2"
	    printf '%s\n' $prefix2 $emtype | ${execdir}/cubeemiss > /dev/null 
	fi
	if [ ! -f ${prefix2}"map-xp-"${emtype}".fits" ]; then
	    echo "Constructing $emtype emissivity maps for $prefix2"
	    printf '%s\n' $prefix2 $emtype | ${execdir}/makemap > /dev/null 
	fi
    fi
done
echo "No more files...done."
exit 0
