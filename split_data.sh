#!/bin/sh

# Get name of file to be read. The master is here
#    /project/data/VariantP/variant_file/homo_sapiens_variation.txt
# but we may want to experiment with a subset of data

if [ $# -lt 1 ]; then
   echo "Usage: readfile filename"
   exit 1
fi
input_file=$1

# Create a top-level directory, under which we will create
# subdirs A,B,C,....

master_dir="/project/scratch/gew123/split_data"

# Count lines as we want to skip header and tail end
# Data is in lines 165-39289586 

nlines=0 

# Bash usually uses spaces to delimit fields, but we want
# entire line read without splitting. Redefine Input Field
# Seperator

OLDIFS=$IFS
IFS=''

while read line
do

  nlines=$((nlines+1))

  if [ $nlines -gt 164 ] && [ $nlines -le 39289586 ]
  then
    # Get Accession code & first character
    AC_Code=`echo $line | awk '{print $2}'`
    first_char=`echo $AC_Code | cut -c1-1`

    # If directory does not exist, create it
    if [ ! -d ${master_dir}/${first_char} ]
    then
      mkdir -p ${master_dir}/${first_char}
    fi

    # Append data into file
    echo $line >> ${master_dir}/${first_char}/${AC_Code}".txt"
  fi

done < ${input_file}

# Return the oroginal Input Filed Seperator (not really
# necessary here, but just good habit)

IFS=$OLDIFS

echo "Lines processed = " $nlines


