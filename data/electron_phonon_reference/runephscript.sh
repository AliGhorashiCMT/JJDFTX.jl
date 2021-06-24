for esmearing in .001 .01 .1 1; do 
	for acceptrange in .001 .01 .1 1; do
		echo esmearing is $esmearing
		echo acceptrange is $acceptrange
		export esmearing
		export acceptrange
		./ephscript.py > ephscript_esmearing"$esmearing"_acceptrange"$acceptrange".out
	done
done
