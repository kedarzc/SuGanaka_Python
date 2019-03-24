import csv
import numpy as np
import math

def readFile(fileName):

	f = open(fileName)
	all_lines = f.readlines()
	f.close()

	return all_lines

def count_occurence_keyword(all_keywords,keyword_to_count):

	# This function counts the occurence of a keyword in the keywords info
	num_keyword_occurence_ids = [i for i, e in enumerate(all_keywords) if e == keyword_to_count]

        # Count the number of occurences
        num_keyword_occurence = len(num_keyword_occurence_ids)
        
	return num_keyword_occurence, num_keyword_occurence_ids


# This function reads the nodal coordinates
def readNodes(all_lines,all_keywords,keyword_lines,all_asterix):

	nativeKeyword = '*Node'

	# How many times does the keyword occur ?
	[numOccurence,k_ids] = count_occurence_keyword(all_keywords,nativeKeyword)

	nodal_keyword_line = keyword_lines[k_ids[0]]
	next_line = all_asterix[all_asterix > nodal_keyword_line].min()

	starting_line = nodal_keyword_line
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
def readElements(all_lines,all_keywords,keyword_lines,all_asterix):

	nativeKeyword = '*Element'

    # How many times does the keyword occur ?
	[numOccurence,k_ids] = count_occurence_keyword(all_keywords,nativeKeyword)

	element_keyword_line = keyword_lines[k_ids[0]]
	next_line = all_asterix[all_asterix > element_keyword_line].min()

	starting_line = element_keyword_line
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

	# This list contains all the keywords in the file
	all_keywords = []
	
	# This list contains all the commented lines 
	comment_lines = []

	# This list contains the line numbers of the keywords
	keyword_lines = []

	with open(fileName) as myFile:

		for num, line in enumerate(myFile, 1):

			if lookup in line and line[1].isalpha() == True:

				current_keyword = line.split(',')[0].rstrip()

				# Save the keywords
				all_keywords.append(current_keyword)

				# List to store the keywords alongwith the 
				# keyword line numbers
				keyword_lines.append(num)

			if lookup in line and line[1] == '*':
				comment_lines.append(num)
        
	# Create a list that has all the lines which start with asterix
	all_asterix = np.array(sorted(keyword_lines + comment_lines))

	return keyword_lines, all_keywords,comment_lines, all_asterix


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

        all_node_nos = Nodes[:,0]
        
	coordinates = np.zeros(shape=(nne,nodof))
	g = np.zeros(shape=(1,nne*nodof))

	# Extract the coordinates of the element
	l = 0
	for k in range(0,nne):
		for j in range(0,nodof):

                        # Node number
                        Node_num_index = np.where(all_node_nos==int(Elements[elemNum][k+1]))
                        
			coordinates[k][j] = Nodes[Node_num_index,j+1]
			g[0,l] = nf[Node_num_index,j]
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

# Create Node Sets
# ----------------

def createNodeSets(all_lines,all_keywords,keyword_lines,all_asterix):

	nativeKeyword = '*Nset'

        # We will create a dictionary to save all  the node sets
        nodeSets = {}
        
	# How many times does the keyword occur ?
	[numOccurence,k_ids] = count_occurence_keyword(all_keywords,nativeKeyword)

        # This loop continues for total number of node sets defined in the input file
        for k in range(0,numOccurence):
                
                # There can be multiple node sets.
                nodal_keyword_line = keyword_lines[k_ids[k]]
                next_line = all_asterix[all_asterix > nodal_keyword_line].min()
                
                starting_line = nodal_keyword_line
                ending_line = next_line

                # Determine the name of the set
                name_of_set =  all_lines[starting_line-1].split('=')[-1]

                my_lines = []
                nodeNums = []

                # Collect the lines that have the nodal numbers
                for i in range(starting_line,ending_line-1):
                         my_lines.append(all_lines[i].split(','))

                # collect the nodal numbers from various lines into ONE list
                for i in range(0,ending_line-1-starting_line):
                        for j in range(0,len(my_lines[i])):
                                nodeNums.append(int((my_lines[i][j]).rstrip()))
                
                nodeSets[name_of_set.rstrip()]= nodeNums
        
	return nodeSets

# This funciton reads the boundary conditions from the input files
def read_BCS(all_lines,all_keywords,keyword_lines,all_asterix):
        
	nativeKeyword = '*Boundary'

	# We will create a dictionary to save all  the node sets
	nodeSets_BCS = {}

	# How many times does the keyword occur ?
	[numOccurence,k_ids] = count_occurence_keyword(all_keywords,nativeKeyword)

	nodal_keyword_line = keyword_lines[k_ids[0]]
	next_line = all_asterix[all_asterix > nodal_keyword_line].min()

	starting_line = nodal_keyword_line
	ending_line = next_line

	# Empty dictionary to save the node sets on which the BCS have been applied
	BCS_NodeSet = {}

	# This set saves the name of the sets on which the boundary condition has been applied
	BC_NodeSet_list = []

	# This set saves the constrained degrees of freedom
	constrained_DOF_string = []

	for i in range(starting_line,ending_line-1):
			BC_NodeSet_list.append((all_lines[i].split(',')[0]).rstrip())
			constrained_DOF_string.append(all_lines[i].split(',')[1::])


	for i in range(0,len(constrained_DOF_string)):
			
		temp_nodal_nums = []

		for j in range(0,len(constrained_DOF_string[i])):
				 temp_nodal_nums.append(int((constrained_DOF_string[i][j]).rstrip()))

		BCS_NodeSet[BC_NodeSet_list[i]] = temp_nodal_nums        

	return BCS_NodeSet


# Apply boundary conditions to the node sets
def apply_BCS(nnd,nodof,Nodes,NodeSets,BCS_NodeSet):

        all_node_nos = Nodes[:,0]
        
        # Initialise the nodal freedom matrix to 1
        nf = np.ones(shape=(nnd,nodof))

        # Determine the total number of BCS
        total_BCS = len(BCS_NodeSet.keys())
        BCS_keys = BCS_NodeSet.keys()

        for i in range(0,total_BCS):
                
                for j in range(0,len(NodeSets[BCS_keys[i]])):

                        # The node number on which the BCS is being applied
                        Node_num = NodeSets[BCS_keys[i]][j]

                        # find the index of the node number
                        Node_num_index = np.where(all_node_nos==Node_num)
                        
                        for k in range(0,len(BCS_NodeSet[BCS_keys[i]])):
                                prescribed_dof = BCS_NodeSet[BCS_keys[i]][k]
                                nf[Node_num_index,prescribed_dof-1] = 0

        return nf
         
        
# Read the loading information
def read_Cloads(all_lines,all_keywords,keyword_lines,all_asterix):

	# This function reads the loads.
	# The type of load is concentrated load on nodes

	nativeKeyword = '*Cload'

	# We will create a dictionary to save all  the node sets
	nodeSets_Cload = {}

	# How many times does the keyword occur ?
	[numOccurence,k_ids] = count_occurence_keyword(all_keywords,nativeKeyword)

	nodal_keyword_line = keyword_lines[k_ids[0]]
	next_line = all_asterix[all_asterix > nodal_keyword_line].min()

	starting_line = nodal_keyword_line
	ending_line = next_line

	# Empty dictionary to save the node sets on which the loads have been applied
	Cload_NodeSet = {}

	# This set saves the name of the sets on which the boundary condition has been applied
	Cload_NodeSet_list = []

	# This set saves the constrained degrees of freedom
	Cload_DOF_MAG_string = []

	for i in range(starting_line,ending_line-1):
		Cload_NodeSet_list.append((all_lines[i].split(',')[0]).rstrip())
		Cload_DOF_MAG_string.append(all_lines[i].split(',')[1::])

	nodal_nums = []

	for i in range(0,len(Cload_NodeSet_list)):

		temp_nodal_nums = []

		for j in range(0,len(Cload_DOF_MAG_string[i])):
			temp_nodal_nums.append(float((Cload_DOF_MAG_string[i][j]).rstrip()))

		nodal_nums.append(temp_nodal_nums)

	return Cload_NodeSet_list, nodal_nums


# This function applied the nodal loads
def apply_cloads(nnd,nodof,Nodes,NodeSets,Cload_NodeSet_list,Cload_dof_mag):

        # Find the node numbers
        all_node_nos = Nodes[:,0]
        
	Nodal_loads = np.zeros(shape=(nnd,nodof))

	# i = node list. Do this loop for all the node lists
	for i in range(0,len(Cload_NodeSet_list)):

	    # node number, do this for all the nodes in the list
	    for j in range(0,len(NodeSets[Cload_NodeSet_list[i]])):

	        # First index tells us the node number in the global nodal_loads array
                # --------------------------------------------------------------------
                # find the index of the node number
                Node_num_index_I = np.where(all_node_nos==NodeSets[Cload_NodeSet_list[i]][j])
	        index_I = Node_num_index_I

	        # Second index tells us the dof of the node in the global nodal_loads array
	        index_II = int(Cload_dof_mag[i][0]-1)

	        Nodal_loads[index_I,index_II] = Cload_dof_mag[i][1]

	return Nodal_loads    

# This function read the material properties from the input file
def read_elastic_material(all_lines,all_keywords,keyword_lines,all_asterix):

	nativeKeyword = '*Elastic'

	# How many times does the keyword occur ?
	[numOccurence,k_ids] = count_occurence_keyword(all_keywords,nativeKeyword)

	nodal_keyword_line = keyword_lines[k_ids[0]]
	next_line = all_asterix[all_asterix > nodal_keyword_line].min()

	starting_line = nodal_keyword_line
	ending_line = next_line

	# Now that we know the starting and ending lines, we will actually read 
	# the elastic modulus and the poisson's ratio
	# We know that the first line will contain the keyword *Elastic and thus 
	# the second line will contain the information we need in the format 
	# E,nu

	for i in range(starting_line,ending_line-1):
		E =  float((all_lines[i].split(',')[0]).rstrip())
		nu =  float((all_lines[i].split(',')[1]).rstrip())

	return E,nu



