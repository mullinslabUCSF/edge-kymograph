function [aftname,forename] = get_image_names(im_dir)
% get_image_names - goes to a directory; identifies all the .tif images;
% and parses the name into two parts: a single forename and a cell string
% array containing all the aftnames. 
% RDM - 7/1/2016

% set the return directory to the current one
ret_dir = pwd;
% go to the specified image directory
cd(im_dir)
% do a linux 'ls' command and store the results in 'pics'. This
% makes pics a list of strings corresponding to the file names.
pics = ls;
% from the list of all files in the directory identify the .tif files and 
% store their positions in 'ind'. Note: in some operating systems (or
% versions of Matlab) the 'ls' and 'dir' commands generate a 2D matrix
% while in others (notably Mac OSX/Matlab 2015b) it generates a 1D vector.
% The following code works only for the latter case.
ind = findstr('.tif',pics);     % find the .tif file position
% find all the line breaks in the string
NLind = findstr(char(10),pics); % ascii code for newline is 10

% make sure pics is a single row vector
pics=pics';
pics=reshape(pics,1,[]);

% use pics to create a list of name parts that can be used to read the 
% original image files or rearranged to construct new names that are linked
% to the originals. This list of name parts gets passed on to the other
% functions
aftname = cell(1,length(ind));
% create a string containing all but the last three digits of the filenames
startname = 1;
for k=1:ind(1)-1    
    % find the first character of the first .tif file name 
    if pics(ind(1)-k)==char(10) || pics(ind(1)-k)==char(9)
        startname = k;
        k = ind(1)-1;
    end
end
forename = pics(1,startname:ind(1)-4);
    
    
for i = 1:length(ind);
    % Convert each element of pics to a simpler name that is just the last
    % three digits of the original name (usually three numbers) plus the
    % '.tif' string. The andor naming convention inserts three zeros, two
    % underscores and the word 'Sequence' into the file name.
    % name{i} = pics(1,ind(i)-20:ind(i)+3); % For Andor produced images
    aftname{i} = pics(1,ind(i)-3:ind(i)+3); % For images cropped by me
end
% in case the files are out of order, sort them by filename
aftname = sort(aftname);
% create a new cell containing all the qdot centroids. This is an array of
% points for every image

cd(ret_dir)

return
end

