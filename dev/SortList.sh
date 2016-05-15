#!/bin/sh
for i in $* 
do
	echo $i	
	new=$i".temp"
	sort +2n -3 +3n -4 $i > $new	
	mv $new $i
done
