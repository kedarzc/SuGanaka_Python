import math

def write_gmsh_file(nameFile,nnd,Nodes,nodesFinal,node_disp,nel,elemName,Elements):

	###
	# VERY IMPORTANT THINGS TO REMEMBER WHEN YOU CHANGE THE CODE
	# 1. The post processing assumes that there are only two steps in the solution
	### 

	nsteps = 2 # Number of steps

	# Create the gmsh file
	gmshfile=open(nameFile,'w')

	# Gmsh version-number file-type data-size
	gmshfile.write('$MeshFormat\n')
	gmshfile.write('2.0 0 8\n')
	gmshfile.write('$End'+'MeshFormat\n')

	# Gmsh Nodal information
	gmshfile.write('$Nodes\n')
	gmshfile.write(str(nnd) + str('\n'))

	for i in range(0,nnd):

		NodeNum = int(Nodes[i][0])
		XCoord = nodesFinal[i][0]
		YCoord = nodesFinal[i][1]
		ZCoord = nodesFinal[i][2]

		gmshfile.write(str(NodeNum) +' ' +str(XCoord) +' '+ str(YCoord) +' '+ str(ZCoord))
		gmshfile.write('\n')

	gmshfile.write('$End'+'Nodes\n')

	# Gmsh Element information
	gmshfile.write('$Elements\n')
	gmshfile.write(str(nel) + str('\n'))
	
	for i in range(0,nel):

		if elemName == 'C3D8':

			ElemNum = int(Elements[i][0])

			N1 = int(Elements[i][1])
			N2 = int(Elements[i][2])
			N3 = int(Elements[i][3])
			N4 = int(Elements[i][4])
			N5 = int(Elements[i][5])
			N6 = int(Elements[i][6])
			N7 = int(Elements[i][7])
			N8 = int(Elements[i][8])

			gmshfile.write(str(ElemNum) +' '+ '5 2 99 2' +' '+ str(N1) +' '+ str(N2) +' '+ str(N3) +' '+ str(N4) +' '+ str(N5) +' '+ str(N6) +' '+ str(N7) +' '+ str(N8))
			gmshfile.write('\n')

	gmshfile.write('$End'+'Elements\n')
	
	# -------------
	# Gmsh NodeData
	# -------------

	displacement_strings = ["U1","U2","U3","U"]


	# Do this for initial step
	for j in range(0,len(displacement_strings)):

		gmshfile.write('$NodeData\n')

		# Type of output
		gmshfile.write(str(1) + '\n' + '"' + str(displacement_strings[j] + '"' + '\n'))

		# Time
		gmshfile.write(str(1) + '\n' + '0.0\n')

		gmshfile.write(str(3) + '\n' + '0\n' + '1\n' + str(nnd) + '\n')

		for i in range(0,nnd):

			NodeNum = i+1
			U_step0 = 0.0
			gmshfile.write(str(NodeNum) + ' ' + str(U_step0) + '\n')

		gmshfile.write('$EndNodeData\n')

	# Repeat for next steps
	for j in range(0,len(displacement_strings)):

		gmshfile.write('$NodeData\n')

		# Type of output
		gmshfile.write(str(1) + '\n' + '"' + str(displacement_strings[j] + '"' + '\n'))

		# Time
		gmshfile.write(str(1) + '\n' + '1.0\n')

		gmshfile.write(str(3) + '\n' + '1\n' + '1\n' + str(nnd) + '\n')

		for i in range(0,nnd):

			NodeNum = i+1

			if displacement_strings[j] == "U1":
				U_stepN = node_disp[i][0]
			
			elif displacement_strings[j] == "U2":
				U_stepN = node_disp[i][1]
			
			elif displacement_strings[j] == "U3":
				U_stepN = node_disp[i][2]
			
			elif displacement_strings[j] == "U":
				U_stepN = math.sqrt((node_disp[i][0]**2 + (node_disp[i][1])**2 + (node_disp[i][2])**2)) 

			gmshfile.write(str(NodeNum) + ' ' + str(U_stepN) + '\n')

		gmshfile.write('$EndNodeData\n')

	# -------------
	# Deformed Mesh
	# -------------

	for n in range(0,nsteps):

		gmshfile.write('$NodeData\n')

		# Type of output
		gmshfile.write(str(1) + '\n' + '"Deformed_Mesh"' + '\n')

		# Time
		gmshfile.write(str(1) + '\n' + str(float(n)) + '\n')

		gmshfile.write(str(3) + '\n' + str(int(n)) +'\n' + '3\n' + str(nnd) + '\n')

		for i in range(0,nnd):

			NodeNum = i+1

			if n == 0:
				gmshfile.write(str(NodeNum) + ' ' + '0.0' + ' ' + '0.0' + ' ' + '0.0' + '\n')
			else:
				gmshfile.write(str(NodeNum) + ' ' + str(node_disp[i][0]) + ' ' + str(node_disp[i][1]) + ' ' + str(node_disp[i][2]) + '\n')

		gmshfile.write('$EndNodeData\n')


		
	gmshfile.close()
