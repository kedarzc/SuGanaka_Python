
def write_gmsh_file(nameFile,nnd,Nodes,nodesFinal,node_disp,nel,elemName,Elements):

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

		gmshfile.write(str(ElemNum) +' '+ '5 2 99 2' +' '+ str(N1) +' '+ str(N2) +' '+ str(N3) +' '+ str(N4) +' '+ str(N1) +' '+ str(N2) +' '+ str(N3) +' '+ str(N4))
		gmshfile.write('\n')

	# gmshfile.write('$End'+'Elements\n')
	# for i in range(0,nel):

	# 	ElemNum = int(Elements[i][0])
	# 	N1 = int(Elements[i][1])
	# 	N2 = int(Elements[i][2])
	# 	N3 = int(Elements[i][3])
	# 	N4 = int(Elements[i][4])

	# 	gmshfile.write(str(ElemNum) +' '+ '3 2 99 2' +' '+ str(N1) +' '+ str(N2) +' '+ str(N3) +' '+ str(N4))
	# 	gmshfile.write('\n')

	# gmshfile.write('$End'+'Elements\n')
	
	# # -------------
	# # Gmsh NodeData
	# # -------------
	# gmshfile.write('$NodeData\n')

	# # Type of output
	# gmshfile.write(str(1) + '\n' + '"U2"\n')

	# # Time
	# gmshfile.write(str(1) + '\n' + '0.0\n')

	# gmshfile.write(str(3) + '\n' + '0\n' + '1\n' + str(nnd) + '\n')

	# for i in range(0,nnd):

	# 	NodeNum = i+1
		
	# 	U1 = node_disp[i][0]
	# 	U2 = node_disp[i][1]
	# 	U3 = 0.0

	# 	gmshfile.write(str(NodeNum) + ' ' + str(U2) + '\n')

	# gmshfile.write('$EndNodeData\n')
	
	gmshfile.close()
