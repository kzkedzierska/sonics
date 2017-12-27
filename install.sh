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
   
   Checks for requirements for SONiCS." 
}

PATH_TO_PYTHON3=$1
if [ "X$PATH_TO_PYTHON3" ==  "X" ]; then
	PATH_TO_PYTHON3=$(which python3);
fi

#check Python version
major=$("$PATH_TO_PYTHON3" --version | cut -f2 -d' ' | cut -f1 -d'.')
minor=$("$PATH_TO_PYTHON3" --version | cut -f2 -d' ' | cut -f2 -d'.')
if [[ major -lt "3" || minor -lt "4" ]]; then
	error "SONiCS requires Python version 3.4 or higher."
fi

#check for all the packages
for package in $(echo "cython logging numpy pandas scipy sklearn pymc argparse shutil"); do
	if ! "$PATH_TO_PYTHON3" -c "import $package" &> /dev/null; then
		missing_packages=$(echo ${missing_packages} $package | sed 's/sklearn/scikit-learn/g');
	fi
done

if [ "X$missing_packages" != "X" ]; then
	error "Missing module(s): ${missing_packages}."
fi

#Cythonize SONiCS
"$PATH_TO_PYTHON3" setup.py build_ext --inplace

#run test 
echo "Running test on support read out: 5|1;6|1;7|5;8|11;9|20;10|24;11|2;12|1"
"$PATH_TO_PYTHON3" run_sonics.py --pvalue_threshold 0.001 --half_random -r 1000 12 "5|1;6|1;7|5;8|11;9|20;10|24;11|2;12|1"

if [ $? == 0 ]; then
    echo "Everythings seems to be working allright "
    echo "The output of the test is in sonics_out.txt"
    rm -rf build
else
    echo "The installation did not work as expected. Please check the error"
    echo "messages above to resolve installation issues. Please contact the" 
    echo "author if you need further help. "
fi
