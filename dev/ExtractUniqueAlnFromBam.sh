samtools view -h $1 | awk '$1 ~/@/ || $5 ==1' | samtools view -bS - > $1.temp
mv $1.temp $2
