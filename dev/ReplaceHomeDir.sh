for i in $*
do
	temp=$i".temp"
	sed 's/\/ifs\/h\/toung/\/home\/jmtoung\/Lab/g' $i > $temp
	mv $temp $i
done 
