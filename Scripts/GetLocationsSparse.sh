
fres=ResultsAll.txt; rm -f $fres
for type in B_1000_Ei R_1000_Ei L_1000_Ei B_1000_1D R_1000_1D L_1000_1D; do
	for rdir in ` ls -d results_SAMC_${type}_Sparse*`; do

   ### location range from MC (PosteriorD.txt)
   fMC=${rdir}/PosteriorD.txt
   if [ ! -e $fMC ]; then
      echo "   Warning: file "$fMC" not found"
      continue
   fi

   # best fitting model
   # check if the last line of fMC is marked as 'best'
	line=`tail -n1 $fMC`
   if [ `echo $line | grep -c 'best'` != 1 ]; then
      echo "best fitting model not found in "$fMC
      continue
   fi
#minfo = (253.971 33.848 -49.665 6.191 245.0761 41.1340 -0.8602) N
	echo $line | awk -F'minfo = \\(' '{print $2}' | awk -v type=$type -F')' '{print $1,type}' >> $fres
	done
	echo -e "\n" >> $fres
done

echo $fres
