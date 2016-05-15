#$ -cwd ### use current directory
#$ -S /bin/bash ### program to execute script
#$ -M toung@mail.med.upenn.edu ### email address
##$ -m ea ### mail is to be sent at abort and end time
#$ -j y ### combine stdout and stderr
#$ -t 1-1 ### array job #
#$ -V ### use current environment variables
##$ -pe DJ 12 ### parallel threads
##$ -l mem_free=10G ### request memory

