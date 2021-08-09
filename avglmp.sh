awk -v k=0 -v u=0 -v e=0 -v p=0 -v n=1 '
{
    k+=$2;u+=$3;e+=$4;p+=$5;printf "%i %.5f %.5f %.5f %.5f\n", n, k/n, u/n, e/n, p/n; n+=1 
}' energieslmp.dat > avglmp.dat
