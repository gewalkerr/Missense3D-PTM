#!/usr/bin/env python

from datetime import datetime

'''
Create directories where the split data will exist
Typically
   cd result_dir
   for i in {0..9};do mkdir $i; done
   for i in {A..Z};do mkdir $i; done
And then in this script set variable result_dir
'''

# Data starts after line number 'start_data' & 
# ends before 'end_data'

data_start = 164
data_end   = 39289587

# Count lines, files, written

line_number = 0
lines_written = 0
files_written = 0
count_tmp = 0

# Input file

input_file = '/project/data/VariantP/variant_data/homo_sapiens_variation.txt'

# Results go here (e.g., esult_dir/A/A0123XYZ.txt). On the linux
# we have already created dir result_dir and sub-dirs 0-9 + A-Z

result_dir = '/project/data/VariantP/split_data'

# Start processing

for line in open(input_file, 'r'):

	line_number += 1
	if (line_number > data_start) and (line_number < data_end):

		words = line.split()
		ac    = words[1]
		lines_written += 1	

		# Just keep track of progress
		count_tmp += 1
		if count_tmp == 10000:
			now = datetime.now()
			time_string = now.strftime("%d/%m/%Y %H:%M:%S")
			print(time_string," : Written ",lines_written," lines")
			count_tmp = 0

		# If first line, record the AC key & open file
		if lines_written == 1:
			ac_old = ac
			filepath = result_dir+'/'+ac[0]+'/'+ac+'.txt'
			files_written += 1
			file_new = open(filepath,'w')

		# Check if AC key changed. If yes, then close old file and open new
		if ac_old != ac:
			file_new.close()
			ac_old = ac
			filepath = result_dir+'/'+ac[0]+'/'+ac+'.txt'
			files_written += 1
			file_new = open(filepath,'w')

		#print(files_written, line_number, ac)
		file_new.write(line)

file_new.close()


