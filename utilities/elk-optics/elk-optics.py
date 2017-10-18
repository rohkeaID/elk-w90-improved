#!/usr/bin/python
#
# Copyright (C) 2015 by Markus Meinert
# This file is distributed under the terms of the GNU General Public License.
# See the file COPYING for license details.
#
# http://elk.sourceforge.net/
#
# elk-optics.py v1.0
#
# This script comes with no warranty. Check the results
# carefully for production use.
#
# Description:
# Reads all EPSILON_ii.OUT files in the present directory
# and computes optical quantities. Diagonal components only.
#
# Input:  EPSILON_ii.OUT 
#
# Output: energy (eV), Re(eps), Im(eps), 
#         refractive index Re(n), Im(n),
#         normal incidence reflectivity R,
#         absorption coefficient alpha (m^-1),
#         EELS -1/Im(eps)
#
# Output is written to optics_ii.out
#

import sys, os, math, cmath

# check which files of type EPSILON_ii.OUT exist and
# return a list of present components
def get_components():
	possible = ['11', '22', '33']
	present = []	
	for p in possible:
		testfilename = 'EPSILON_%s.OUT' % p
		if os.path.isfile(testfilename):
			present.append(p)
	return present

# read the EPSILON_ii.OUT file
# return lists of energies and complex eps
def read_epsilon(filename):
	handle = open(filename, 'r')
	content = handle.readlines()
	handle.close()
	
	data = [[],[]]
	for line in content:
		l = line.split()
		if l == []:
			continue
		data[0].append(float(l[0]))
		data[1].append(float(l[1]))

	# energies are read from first column of the first data block
	# real part of epsilon is read from the second column of the first data block (first half of the data)
	# imaginary part of epsilon is read from the second column of the second data block (second half of the data)
	datalength = int( len( data[0] ) / 2.)
	energies = data[0][0:datalength]
	eps_cplx = [complex(a,b) for a,b in zip(data[1][0:datalength], data[1][datalength:])]

	return energies, eps_cplx

# compute optical properties from energies and complex epsilon
def write_optical_properties(energies, eps_cplx, component):
	# complex refractive index N and extinction coefficient kappa
	# complex refractive index: N = n_r + ik
	N = [cmath.sqrt(x1) for x1 in eps_cplx]
	k = [cmath.sqrt(x1).imag for x1 in eps_cplx]

	# normal incidence reflectivity from complex refractive index
	R = [abs((1.-x1)/(1.+x1))**2 for x1 in N]

	# absorption coefficient in SI units from extinction coefficient and energy
	Ha_to_J = 27.21138602 * 1.6021766208E-19
	hbar = 6.626070040E-34 / (2 * math.pi)
	c = 2.99792458E8
	Ha_to_omegaSI = Ha_to_J / hbar

	alpha = [2 * (x1 * Ha_to_omegaSI) / c * x2 for x1, x2 in zip(energies, k)]
	
	# format data and write to file optics_ii.out
	data = zip(energies, eps_cplx, N, R, alpha)
	output = '%14s %14s %14s %14s %14s %14s %14s %14s\n' % ('#  energy (eV)', 'Re(eps)', 'Im(eps)', 'Re(n)', 'Im(n)', 'R', 'alpha (m^-1)', 'EELS')
	for line in data:
		output += '%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n' % (line[0]*27.21138602, line[1].real, line[1].imag, line[2].real, line[2].imag, line[3], line[4], -(1/line[1]).imag)
	outfilename = 'optics_%s.out' % component
	outfile = open(outfilename, 'w')
	outfile.write(output)
	outfile.close()

# main loop over diagonal components of the epsilon tensor
print('===================')
print('| elk-optics v1.0 |')
print('===================')
print
print('Looking for EPSILON_ii.OUT files...')

components = get_components()

if components == []:
	sys.exit('No EPSILON_ii.OUT files found. Exit.\n')
else:
	print('Files found:')
	for c in components:
		print('    EPSILON_%s.OUT') % c
	print

for c in components:
	filename = 'EPSILON_%s.OUT' % c
	print('Working on %s ...') % filename
	energies, eps_cplx = read_epsilon(filename)
	write_optical_properties(energies, eps_cplx, c)
	print('Optical properties written to optics_%s.out') % c
	print
