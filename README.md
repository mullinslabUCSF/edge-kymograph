# edge-kymograph
Matlab code for tracking edge of cell and extracting kymographs

Droid5

About
A function in Matlab that is designed to find the moving perimeter of a cell and follow protein dynamics around it. The input is the name of the directory containing a numbered sequence of .tif images of up to 2 marker proteins at the cell edge. The output is a space-time kymograph displaying the dynamics and trajectories of the edge protein(s) of interest.

The following accessory functions in addition to droid5.m are required:
get_image_names.m
bdy_issues.m
boundary_lines.m
generic_edge2.m
makeline.m
pathfinder2.m

Compatibility
Written for Matlab R2019a.

Example
We originally developed Droid5 to analyze the initiation and disassembly dynamics of leading-edge clusters containing the proteins VASP and lamellipodin (Lpd) in B16F1 cells. See preprint for more details. 

Example input below is a two-color fluorescent image of leading-edge proteins VASP (left, green) and Lpd (right, magenta) from a time-lapse movie. 
