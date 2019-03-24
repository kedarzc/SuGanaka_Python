import StandardLibrary as STDLIB
import preProcess as PREPRO


# INPUTS 
# -----------------------------------------------------------------------------------------
# The name of the input file
input_file = 'Beam_bending.sinp'


# -----------------------------------------------------------------------------------------
# INPUTS END

name_ip_file = input_file[:-5]

# Read entire input file
all_lines = STDLIB.readFile(input_file)

# Parse the input file for 'keywords'
[keyword_lines, all_keywords,comment_lines, all_asterix] = STDLIB.parseKeywords(input_file)

# Read the Nodes
Nodes = STDLIB.readNodes(all_lines,all_keywords,keyword_lines,all_asterix)

# Read the elements
Elements = STDLIB.readElements(all_lines,all_keywords,keyword_lines,all_asterix)

# Count the number of nodes
nnd = len(Nodes)
nel = len(Elements)

PREPRO.write_mesh(name_ip_file,nnd,nel,Nodes,Elements)