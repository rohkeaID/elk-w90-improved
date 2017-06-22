#!/usr/bin/python

# elk blocks to columns by Markus Meinert
# update January 2015, fixed a bug related to blank lines

# Usage: blocks2columns.py PDOS_S01_A0001.OUT
#        blocks2columns.py TDOS.OUT
#        blocks2columns.py BAND.OUT

import sys, os

print("\n =========================\n | elk blocks to columns |\n =========================\n")

# Read the file.
filename = sys.argv[1]
f = open(filename, 'r')
data = f.readlines()
f.close()

# Analyze the file.

# Number of lines.
nlines = len(data)
print(" Number of lines: %i " % nlines)

# Count blank lines to determine number of datasets.
ndatasets = 0
for line in data:
	if line.split() == []:
		ndatasets += 1
# If last line is not blank, add one to ndatasets and nlines.
if data[-1].split() != []:
	ndatasets += 1
	nlines += 1
print(" Number of datasets: %i " % ndatasets)

# Number of lines per block is:
nlinesperblock = (nlines - ndatasets)/ndatasets
print(" Number of lines per block: %i " % nlinesperblock)

# Collect the datasets into a list of lists with a double-loop over datasets and lines.
datasets = []
for i in range(0,ndatasets):
        currentset = []
        for j in range(i*nlinesperblock + i, (i*nlinesperblock + i) + nlinesperblock):
                currentset.append(data[j].split()) # Split each line by empty spaces.
        datasets.append(currentset)

output = ""

# Generate a head line
output += "#%21s" % "x-axis"
for i in range(1,ndatasets+1):
	blockname = "block_%i" % i
	output += "%22s" % blockname
output += "\n"

# Merge the datasets line-wise.
for i in range(0,nlinesperblock):
	# x-axis as first column, read from first block
	line = '%22.13e' % (float(datasets[0][i][0]))
	# Append the block values as columns.
	for j in range(0, ndatasets):
		line += '%22.13e' % (float(datasets[j][i][1]))
	line += "\n"
	output += line

filename = filename + ".columns"
if os.path.exists(filename):
	print("\n Output file %s exists. Exit.\n" % filename)
else:
	f = open(filename, 'w')
	f.write(output)
	f.close()
	print("\n Output filename: %s\n Done.\n" % filename)
