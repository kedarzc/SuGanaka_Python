import csv
import numpy as np
import math

def readFile(fileName):

	f = open(fileName)
	all_lines = f.readlines()
	f.close()

	return all_lines


# This function reads the nodal coordinates
def readNodes(all_lines,keyword_info,all_keywords_line_nos,all_asterix):

	# print keyword_info

	nodal_keyword = keyword_info['*Node']
	next_line = all_asterix[all_asterix > nodal_keyword].min()

	# print next_line

	starting_line = nodal_keyword
	ending_line = next_line
		
	Nodes = np.zeros(shape=(ending_line-starting_line-1,3))

	Nodal_counter = 0

	for i in range(starting_line+1,ending_line):
		Nodes[Nodal_counter,0] = float(all_lines[i-1].split(',')[0])
		Nodes[Nodal_counter,1] = float(all_lines[i-1].split(',')[1])
		Nodes[Nodal_counter,2] = float(all_lines[i-1].split(',')[2])
		
		Nodal_counter = Nodal_counter + 1

	return Nodes



# This function reads the nodal coordinates
def readElements(all_lines,keyword_info,all_keywords_line_nos,all_asterix):

	# print keyword_info

	nodal_keyword = keyword_info['*Element']
	next_line = all_asterix[all_asterix > nodal_keyword].min()

	# print next_line

	starting_line = nodal_keyword
	ending_line = next_line
		
	Elements = np.zeros(shape=(ending_line-starting_line-1,5))

	Element_counter = 0

	for i in range(starting_line+1,ending_line):
		Elements[Element_counter,0] = float(all_lines[i-1].split(',')[0])
		Elements[Element_counter,1] = float(all_lines[i-1].split(',')[1])
		Elements[Element_counter,2] = float(all_lines[i-1].split(',')[2])
		Elements[Element_counter,3] = float(all_lines[i-1].split(',')[3])
		Elements[Element_counter,4] = float(all_lines[i-1].split(',')[4])
		
		Element_counter = Element_counter + 1

	return Elements




def parseKeywords(fileName):

	lookup = '*'

	all_keywords = []
	all_keywords_line_nos = []
	comment_lines = []


	keyword_info = {}

	with open(fileName) as myFile:
		for num, line in enumerate(myFile, 1):

			if lookup in line and line[1].isalpha() == True:

				current_kewyword = line.split(',')[0]

				# Dictionary to store the keywords alongwith the 
				# keyword line numbers
				keyword_info[current_kewyword.rstrip()] = num

			if lookup in line and line[1] == '*':
				comment_lines.append(num)

	# Create array from the line numbers
	all_keywords_line_nos = np.array(keyword_info.values())

	# Create array from the keywords
	all_keywords = keyword_info.keys()

	# Create an a list that has all the lines which start with asterix
	all_asterix = np.array(sorted(keyword_info.values() + comment_lines))

	return keyword_info, all_keywords, all_keywords_line_nos,comment_lines, all_asterix


def gaussPoints(ng):

	# This function returns the abcissae and the weights of the Gauss points 
	# for number of gauss points

	samp = np.zeros(shape=(ng,2))

	# Nomenclature --> samp[abcissa,weight]

	if ng == 1:

		samp[0,0] = 0.0

		samp[0,1] = 2.0

	if ng == 2:

		samp[0,0] = -1.0/math.sqrt(3.0) 
		samp[1,0] = 1.0/math.sqrt(3.0)

		samp[0,1] = 1.0 
		samp[1,1] = 1.0

	if ng == 3:

		samp[0,0] = -2.0/math.sqrt(15.0)
		samp[1,0] = 0.0
		samp[2,0] = +2.0/math.sqrt(15.0)
 
		samp[0,1] = 5.0/9.0
		samp[1,1] = 8.0/9.0
		samp[2,1] = 5.0/9.0

	return samp


# Extract the coordinates for the Q4 element and 
# the steering vector for the element
def elem_Q4(elemNum,Nodes,Elements,nne,nodof,nf):

	coordinates = np.zeros(shape=(nne,nodof))
	g = np.zeros(shape=(1,nne*nodof))

	# Extract the coordinates of the element
	l = 0
	for k in range(0,nne):
		for j in range(0,nodof):
			coordinates[k][j] = Nodes[int(Elements[elemNum][k+1])-1,j+1]
			g[0,l] = nf[int(Elements[elemNum][k+1])-1,j]
			l = l + 1

	return coordinates, g



# This function returns the shape functions (fun)
# and their derivatives w.r.t. xi and eta at each Gauss Point
def fmQ4_lin(samp,ig,jg):

	xi  = samp[ig][0]
	eta = samp[ig][1]

	# Make this arbitrary
	fun = np.zeros(shape=(4,1))
	der = np.zeros(shape=(2,4))

	# Form the shape functions
	fun[0][0] = 0.25 * (1.0 - xi - eta + xi*eta)
	fun[1][0] = 0.25 * (1.0 + xi - eta - xi*eta)
	fun[2][0] = 0.25 * (1.0 + xi + eta + xi*eta)
	fun[3][0] = 0.25 * (1.0 - xi + eta - xi*eta)

	# Form the derivative of shape functions
	der[0][0] = 0.25 * -(1.0-eta)
	der[0][1] = 0.25 * +(1.0-eta)
	der[0][2] = 0.25 * +(1.0+eta)
	der[0][3] = 0.25 * -(1.0+eta)

	der[1][0] = 0.25 * -(1-xi)
	der[1][1] = 0.25 * -(1+xi)
	der[1][2] = 0.25 * +(1+xi)
	der[1][3] = 0.25 * +(1-xi)

	return der, fun


# This function forms the B-matrix of a linear Q4 element
def formbee_Q4_lin(deriv,nne,eldof):

	bee = np.zeros(shape=(3,eldof))

	# First row
	bee[0][0] = deriv[0][0]
	bee[0][2] = deriv[0][1]
	bee[0][4] = deriv[0][2]
	bee[0][6] = deriv[0][3]

	# Second row
	bee[1][1] = deriv[1][0]
	bee[1][3] = deriv[1][1]
	bee[1][5] = deriv[1][2]
	bee[1][7] = deriv[1][3]

	# Third row
	bee[2][0] = deriv[1][0]
	bee[2][1] = deriv[0][0]
	bee[2][2] = deriv[1][1]
	bee[2][3] = deriv[0][1]
	bee[2][4] = deriv[1][2]
	bee[2][5] = deriv[0][2]
	bee[2][6] = deriv[1][3]
	bee[2][7] = deriv[0][3]

	return bee

# Form the global stiffness matrix
def form_KK(KK,ke,g,eldof):

	for i in range(0,eldof):
		if g[0][i] != 0:
			for j in range(0,eldof):
				if g[0][j] != 0:
					KK[int(g[0][i]-1),int(g[0][j]-1)] = KK[int(g[0][i])-1,int(g[0][j]-1)] + ke[i][j]

	return KK

# Seperate displacements
def seprarate_disp(nodof,nnd,delta,nf):

	node_disp = np.zeros(shape=(nnd,nodof))

	for j in range(0,nodof):
		for i in range(0,nnd):
			if nf[i][j] != 0:
				node_disp[i][j] = delta[int(nf[i][j])-1]

	return node_disp 
