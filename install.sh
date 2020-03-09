#!bin/bash
PROGRAM=$(basename $0)

error()
{
    echo -e "$@" 1>&2
    exit 1
}

usage()
{
   echo -e "Usage:

   $PROGRAM [PATH_TO_PYTHON3]
   
   Checks for requirements and compiles SONiCS. Performs small test to check if everything is working properly." 
}

PATH_TO_PYTHON3=$1
if [ "X$PATH_TO_PYTHON3" ==  "X" ]; then
	PATH_TO_PYTHON3=$(which python3);
fi

#check Python version
major=$("$PATH_TO_PYTHON3" --version 2>&1 | cut -f2 -d' ' | cut -f1 -d'.')
minor=$("$PATH_TO_PYTHON3" --version 2>&1 | cut -f2 -d' ' | cut -f2 -d'.')
if [[ major -lt "3" || minor -lt "4" ]]; then
	error "SONiCS requires Python version 3.4 or higher."
fi

#check for all the packages
for package in $(echo "cython numpy pandas scipy pymc"); do
	if ! "$PATH_TO_PYTHON3" -c "import $package" &> /dev/null; then
		missing_packages=$(echo ${missing_packages} $package);
	fi
done

if [ "X$missing_packages" != "X" ]; then
	error "Missing module(s): ${missing_packages}."
fi

#Cythonize SONiCS
"$PATH_TO_PYTHON3" setup.py build_ext --inplace

#run test 
echo "Running test on support read out: 8|5;9|113;10|89"
"$PATH_TO_PYTHON3" sonics "8|5;9|113;10|89"

if [ $? == 0 ]; then
    echo "Everything seems to be working all right."
    echo "The output of the test is in sonics_out.txt"
    rm -rf build
else
    echo "The installation did not work as expected. Please check the error"
    echo "messages above to resolve installation issues. Please contact the" 
    echo "author if you need further help. "
fi
