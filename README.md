# getMesh
Function to gather mesh information from gmesh files



### DEFINE SHORT-NAME OUTPUTS:

O1=msh; Mesh file from gmesh
O2=Mdpts; Midpoint Database
O3=cons; Various important constants
O4=solidpts; Nodes (vertices and midpoints) on solid domain
O5=fluidpts; Nodes (vertices and midpoints) on fluid domain
O6=solperimeter; Nodes (vertices and midpoints) on solid/fluid bdry
O7=outerbdrynodes; Vertices on the outer boundary
O8=ELsolid; Global element numbers of vertices on solid domain
O9=ELfluid; Global element numbers of vertices on fluid domain
O10=edgetri; List of global element #'s in solid touching bdry
