Initiation and disassembly of filopodia tip complexes containing VASP and lamellipodin
Karen W. Cheng, R. Dyche Mullins
MBoC 2020
https://doi.org/10.1091/mbc.E20-04-0270

# edge-kymograph
Matlab code for tracking edge of cell and extracting kymographs

**Droid5**

**About**

A function in Matlab that is designed to find the moving perimeter of a cell and follow protein dynamics around it. The input is the name of the directory containing a numbered sequence of .tif images of up to 2 marker proteins at the cell edge. The output is a space-time kymograph displaying the dynamics and trajectories of the edge protein(s) of interest.

The following accessory functions in addition to droid5.m are required:

get_image_names.m

bdy_issues.m

boundary_lines.m

generic_edge2.m

makeline.m

pathfinder2.m

**Compatibility**

Written for Matlab R2019a.

**Example**

We originally developed Droid5 to analyze the initiation and disassembly dynamics of leading-edge clusters containing the proteins VASP and lamellipodin (Lpd) in B16F1 cells. See [preprint](https://www.biorxiv.org/content/10.1101/2020.02.21.960229v1) for more details. 

Example input below is a two-color fluorescent image of leading-edge proteins VASP (left, green) and Lpd (right, magenta) from a time-lapse movie. 

![Image 0](https://github.com/mullinslabUCSF/edge-kymograph/blob/master/images/image_0.jpg?raw=true)

Example outputs below are (1) overlaid traces of the cell edge perimeter region of interest color coded by time frame and (2) adaptive kymographs of protein cluster dynamics of  both VASP and Lp channels.

![Image 1](https://github.com/mullinslabUCSF/edge-kymograph/blob/master/images/image_1.jpg?raw=true)

The data set corresponding to this example is available in the folders ‘ch1_VASP’ and ‘ch2_LPD’.

**Workflow**

1. Split the two-color time lapse movie into separate channels and save each as an ‘image sequence’ in separate folders.

2. Create new folders in the same directory called ‘results’ and ‘masked’ to store the outputs files.

![Image 2](https://github.com/mullinslabUCSF/edge-kymograph/blob/master/images/image_2.jpg?raw=true)

3. Set the parameter values (blur, threshold, ch1_maxintensity, ch2_maxintensity) in path_param. These values will differ depending on the fluorophore, levels of protein expression, etc and should be checked for each cell that is analyzed. 

    1. For reference, the path_param values used to analyze the sample data set is path_param = [5, 0.0049, 4000, 4000].  

    2. You can separately run the bdy_issues.m function to determine appropriate threshold and blurring values that yield accurate cell boundary lines before running the main droid5.m function. This function takes in threshold and blur values and outputs a figure with the boundary lines for the image sequence overlaid. 

4. Once good path_param values are established, run the main function: droid5.m.

5. An interactive figure will pop up with the overlaid boundary lines of the cell. To help the function determine the region of interest, click on the inside and outside of the (1) start point (2) end point (3) midpoint of the boundary lines. Then click twice outside of the cell to make a small square (4) for background subtraction for photobleaching calculations. 

![Image 3](https://github.com/mullinslabUCSF/edge-kymograph/blob/master/images/image_3.jpg?raw=true)

6. Once the dynamic region of interest is selected with input from the user, droid5 runs and outputs the following figures and saves the kymograph information to the results directory. 

![Image 4](https://github.com/mullinslabUCSF/edge-kymograph/blob/master/images/image_4.jpg?raw=true)
