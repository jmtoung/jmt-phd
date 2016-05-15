#!/bin/sh
for i in $*
do 
	echo $i 
	new=$i".temp"
	sort +0 -1 +1n -2 +2 -3 +3 -4 +4 -5 +5n -6 +6n -7 +7n -8 $i > $new
	mv $new $i
done
