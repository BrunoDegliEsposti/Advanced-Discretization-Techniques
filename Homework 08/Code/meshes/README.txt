Contains meshes of different 2D geometries:
 * square
 * L-shaped domain
 * disk

For each geometry, meshes with different refinement levels are 
provided, the filenames have the format
  *geometry*-*refinement level*.msh
e.g. "lshaped-2.msh".

You do not need to change anything in this directory.

If you wish to explore other geometries, or refine further, you can 
generate you own .msh files using the original source files (*.geo).
They describe the outer boundaries of the mesh and can be opened 
with the open source mesh generation software GMSH
  http://www.gmsh.info
