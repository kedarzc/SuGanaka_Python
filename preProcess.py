
def write_mesh(nameFile,nnd,nel,Nodes,Elements):

	# Write the mesh file
	# -------------------

	# Create the gmsh file
	gmshfile=open(nameFile + '_mesh.msh','w')

	# Gmsh version-number file-type data-size
	gmshfile.write('$MeshFormat\n')
	gmshfile.write('2.0 0 8\n')
	gmshfile.write('$End'+'MeshFormat\n')

	# Gmsh Nodal information
	gmshfile.write('$Nodes\n')
	gmshfile.write(str(nnd) + str('\n'))

	for i in range(0,nnd):

		NodeNum = i+1
		XCoord = Nodes[i][1]
		YCoord = Nodes[i][2]
		ZCoord = 0.0

		gmshfile.write(str(NodeNum) +' ' +str(XCoord) +' '+ str(YCoord) +' '+ str(ZCoord))
		gmshfile.write('\n')

	gmshfile.write('$End'+'Nodes\n')

	# Gmsh Element information
	gmshfile.write('$Elements\n')
	gmshfile.write(str(nel) + str('\n'))

	for i in range(0,nel):

		ElemNum = int(Elements[i][0])
		N1 = int(Elements[i][1])
		N2 = int(Elements[i][2])
		N3 = int(Elements[i][3])
		N4 = int(Elements[i][4])

		gmshfile.write(str(ElemNum) +' '+ '3 2 99 2' +' '+ str(N1) +' '+ str(N2) +' '+ str(N3) +' '+ str(N4))
		gmshfile.write('\n')

	gmshfile.write('$End'+'Elements\n')

	gmshfile.close()
