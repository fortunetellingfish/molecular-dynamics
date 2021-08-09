for x in "$@" #loop over all input files
do
touch temp.dat

awk '{if(NR==1){ d=$1; n=$2; print d "    " n "    "  0;} else print $0}' $x >> temp.dat #use the file provided in command line argument 1
mv temp.dat $x
done
