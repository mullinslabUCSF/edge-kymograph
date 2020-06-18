function [paths,pimage] = pathfinder2(imgdir,result_dir,name,forename,thresh,blur,si,ei,mi)
% pathfinder() - This function identifies the edge of a cell expressing a
% fluorescent protein. This edge becomes the ROI for identifying particles,
% such as VASP clusters. 
% adapted from 'generic_edge()' by RDM 5/3/2019
%
% paths - is a cell array with leading-edges paths for each cell
%
% imgdir - is the directory that contains the timelapse image sequence  
% name - is a vector of name suffixes for reconstrucitng file names
% forename - is the name prefix for reconstructing file names
% thresh - is the threshold used to discriminate objects from noise
% blur - is the degree of gaussian blurring to perform on initial image
% si,ei,mi - are binary images used to find start, end, and mid points
% bi - is a binary image used to compute the background noise in the image

% set the return directory to the current one
ret_dir = pwd;

% count the number of frames
n_fr = length(name);
% initialize cell arrays to hold the 'path' data
% the first position is the number of points (n)
% the second position is the nx2 matrix of x,y coordinates of the cell edge
paths = cell(n_fr,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up an RGB image to hold the output cell boundary and ROI image
% first open the first image in the stack and use it to get parameters
cd(imgdir)
I = imread([forename,name{1}]);
ipxy = size(I);
% initialize the RGB image
pimage = double(ones(ipxy(1),ipxy(2),3));
% create a color map for drawing the leading edge
pcmap = jet(64);
cd(ret_dir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lwid - is the width of the region of interest at the cell edge
lwid = 20;
lrad = floor(lwid/2);

% run through all the frames in the directory and build the ROI's
for i=1:n_fr
    
    % go to the channel 1 image directory 
    cd(imgdir)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    I = imread([forename,name{i}]);
    % blur the image 
    I2 = imgaussfilt(I,blur);
    %Median Filter
    I3 = medfilt2(I2,[1 1]);
    %Subtract Background
    %I4 = imsubtract(I3,imopen(I3,strel('disk',2))); % last # is the radius of the objects
    I4=I3;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % the bw variables hold sequential, processed versions of the original
    % image. The first step (bw) sets a threshold
    % the second step (bw2) fills in any 'holes' in the image
    % the third step (bw3) determines the size of the blobs kept
    % the fourth step (bw4) creates a (smoothed?) perimeter around every image
    % the fifth step (bw5) fills the perimeters back in
%     bw = im2bw(I4,thresh);
    bw = imbinarize(I4,thresh);
    bw2 = imfill(bw,'holes'); 
    % may not need to use bwareaopen
    % but bw3 and bw4_perim are the two most useful bw images for analysis
    bw3 = bwareaopen(bw2,100); 
    bw4_perim = bwperim(bw3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    hold on
    % dilate the detected perimeter to capture all the leading edge signal
    SE = strel('disk', lrad , 4);
    bw5 = imdilate(bw4_perim,SE);
    
    % use the edge mask (bw5) to mask the original image (I)
    I = I.*uint16(bw5);   
    
    % figure out where to break the perimeter for the kymograph
    tempi = si.*bw4_perim;
    [tmx,tmy] = find(tempi);
    if isempty(tmx)
        % this desperate measure will not work on the first pass, because
        % br_pt has not been set yet. It should work on subsequent passes
        [br_pt(1),br_pt(2),k]=findthread(bw4_perim,br_pt(1),br_pt(2),100);
    else
        br_pt = [tmx(1) tmy(1)];
    end
    % find the endpoint of the ROI
    tempi = ei.*bw4_perim;
    [tmx,tmy] = find(tempi);
    if isempty(tmx)
        % this desperate measure will not work on the first pass, because
        % end_pt has not been set yet. It should work on subsequent passes
        [end_pt(1),end_pt(2),k]=findthread(bw4_perim,end_pt(1),end_pt(2),100);
    else
        end_pt = [tmx(1) tmy(1)];
    end
    % find the midpoint of the ROI
    tempi = mi.*bw4_perim;
    [tmx,tmy] = find(tempi);
    if isempty(tmx)
        % this desperate measure will not work on the first pass, because
        % mid_pt has not been set yet. It should work on subsequent passes
        [mid_pt(1),mid_pt(2),k]=findthread(bw4_perim,mid_pt(1),mid_pt(2),100);
    else
        mid_pt = [tmx(1) tmy(1)];
    end
    
    % on the first pass set up an image for overlaying ROIs
    if i==1
        figure(2);
        JJ = imcomplement(bw4_perim);
        imshow(JJ);
        % also draw the outline of all the cells in the field into pimage
        pimage(:,:,1) = pimage(:,:,1) - bw4_perim;
        pimage(:,:,2) = pimage(:,:,2) - bw4_perim;
        pimage(:,:,3) = pimage(:,:,3) - bw4_perim;
    end    
    
    % order the points around the cell perimeter
    [path1]=walkabout(bw4_perim,br_pt);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % draw the new cell outline (in black) to the output RGB image
    for ip = 1:length(path1)
        pimage(path1(ip,1),path1(ip,2),:)=0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % find the end and mid points and see whether the path is going in the
    % right direction (from br -> mid -> end)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ugly ugly ugly - but works (mostly)   %
    % start with vectors to hold distances  %
    % between end or mid points and path    %
    rmid=0;
    rend=0;
    for k=1:length(path1)                   %
        middist=(path1(k,1)-mid_pt(1))^2 + (path1(k,2)-mid_pt(2))^2;%
        enddist=(path1(k,1)-end_pt(1))^2 + (path1(k,2)-end_pt(2))^2;%
        if k==1
            oldmid = middist;
            oldend = enddist;
        else
            if middist<oldmid
                oldmid = middist;
                rmid=k;
            end
            if enddist<oldend
                oldend = enddist;
                rend=k;
            end
        end
    end                                     %
    
    if rmid>rend
        % flip the path and recalculate rend
        rend=length(path1)-rend; 
        path1=flipud(path1);
    end
    % set the path to be the subset of the perimeter
    pathe=path1(1:rend,:);
    num=length(pathe);
    
    %%%%%%% temporary patch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fix occasional artifact where starting point is repeated
    ptrim = 0;
    for ipt = 2:num
        if pathe(ipt,1)==pathe(ipt-1,1) && pathe(ipt,2)==pathe(ipt-1,2)
            ptrim = ptrim + 1;
        end
    end
    if ptrim > 0
        pathe = pathe(1+ptrim:num,:);
    end
    num = length(pathe);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % load the paths cell array
    paths{i,1}=num;
    paths{i,2}=pathe;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % draw the new ROI outline on the output RGB image in a new color
    cmap=parula(64);
    oin=max([1 floor((i/n_fr)*64)]);
    %pcolor = [(i/n_fr) (i/n_fr) 1-(i/n_fr)];
    for ip = 1:num
        pimage(pathe(ip,1),pathe(ip,2),:)=cmap(oin,:);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % add the new ROI outline to the existing figure 1
    figure(2);
    cmap=parula(64);
    oin=max([1 floor((i/n_fr)*64)]);
    hold on
    plot(path1(:,2),path1(:,1),'color',[0 0 0]);
    plot(pathe(:,2),pathe(:,1),'color',cmap(oin,:));
    hold off
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now save the cell image that was just generated
cd(result_dir)
% grab the image from the figure
F = getframe(gcf);
% save the file with a good filename
filename = ['outline_ROI.jpg'];
imwrite(F.cdata,filename);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% there's no place like home
cd(ret_dir)

return
end

function [path]=walkabout(inbwimage,start)
% this function takes a binary image of a closed (Jordan) curve and
% returns an ordered array of x,y indices that walks, stepwise around the
% entire curve

% path - is a two-column matrix containing the ordered x and y coordinates
% num - is the number of points in the closed curve
% bwimage - is the input image containing the one-pixel-wide curve
% start - is a pair of coordinates defining the starting position
%           N.B. if start = [0 0] the alogrithm automatcally breaks
%           symmetry at a point near the middle left edge of the cell

% maxrad - is the maximum search radius allowable for finding the next
% point in the perimeter. This is passed to the function findthread() and
% it is used to pad the original image. For the moment the value of maxrad
% is HARD CODED.
maxrad = 20;
% first, we build a padded version of the input image that has a boundary
% buffer of zeros whose width is equal to maxrad. This enables the code to 
% better deal with edge crashes. Specifically, we will embed the original 
% nxn image inside a larger,(n+2k)x(n+2k), image
orgsize = length(inbwimage);         % size of input image
padsize = 2*maxrad+orgsize;           % size of padded image
bwimage = zeros(padsize,padsize);   % initialize padded image
% load original image into padded holder
bwimage(maxrad+1:maxrad+orgsize,maxrad+1:maxrad+orgsize)=inbwimage;
% remove unconnected 'spurs'
%%%%%%%%%%%% NEW CODE GOES HERE %%%%%%%%%%%%%%%%%%%%%%%


% first find all of the non-zero elements in the image (row and column are
% the indices of the non-zero elements; v is the value of each pixel)
[row,col,v]=find(bwimage);

% estimate the length of the output matrix
num = length(row);
% at this point we don't know the output matrix length because there might
% be spurious stuff in the image. Currently num is a useful upper bound.

% initialize vectors to hold the output (sorted) rows and columns
rf = zeros(num,1);
cf = zeros(num,1);

% create starting points corrected for the added padding
rp = start(1) + maxrad;
cp = start(2) + maxrad;
% store the starting point in the output vectors
rf(1) = start(1);
cf(1) = start(2);

% now run through the connected points on the curve and order them
i=1;
while i<=num-1
    % remove the point at the current position
    if i==1
        % to prevent hopping back and forth across the initial break point
        % we clear a 12x12 pixel zone around the starting point
        bwimage(rp-12:rp+12,cp-12:cp+12)=zeros(25);
    else
        bwimage(rp,cp) = 0; 
    end
    % use findthread() to get the next point
    [rp,cp,k]=findthread(bwimage,rp,cp,maxrad);
    if rp ~= 0 && cp ~= 0
        % store the new buffer-corrected boundary points in the output
        % vectors
        rf(i+1)=rp-maxrad;
        cf(i+1)=cp-maxrad;
        % increment the index
        i=i+1;
    else
        % set the number of points detected in the boundary
        num = i-1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the output parameters
% take account of the fact that the trace probably ends early
path = [rf(1:num) cf(1:num)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return
end

function [startpt,endpt,midpt]=boundary_set(bwimage)
% Function boundary_set enables the user to define a region of interest 
% along the perimeter of the input image(bwimage) using the mouse. 
% The function returns the coordinates of the start and end points together
% with a midpoint required to define which segment of the perimeter to
% analyze.

% First we invert the image pixel values so that the background is white
% and the cell perimeter is black. This gives better contrast for clicking
% because the cursor is also black.
figure(1);
bwimage2=imcomplement(bwimage);
text_str = 'Click three points to define ROI: (1)start, (2)end, (3)midpoint';
imshow(bwimage2)
text(15,15,text_str);

% initialize the data values
x = zeros(3,1);
y = zeros(3,1);

% set up a colormap for plotting the breakpoints
cmap=hot(8);
hold on
% now get three points from the image
for i=1:3
    [x(i),y(i)]=ginput(1);
    % make sure the values are well behaved
    y(i)=floor(y(i));
    x(i)=floor(x(i));
    % now plot the chosen point on the image
    rectangle('Position',[x(i)-5 y(i)-5 10 10],'Curvature',[1 1],'EdgeColor',cmap(i,:));
end
hold off

% load up the return variables
startpt=[y(1),x(1)];
endpt=[y(2),x(2)];
midpt=[y(3),x(3)];
% because in matlab, [col row] = [y,x]

return
end

function [rf,cf,k]=findthread(bwimage,rp,cp,maxk)
% This function searches the input image 'bwimage' for the nearest nonzero
% element to the input row and column parameters 'rp' and 'cp'. The outputs
% are the image point and size of the search submatrix that identified it
% NOTE THAT THIS CODE ASSUMES THE IMAGE IS SQUARE

% rp and cp - are the input row and coumn position for the search
% bwimage - is a logical image containing the perimeter of an object
% maxk - is the maximum search radius allowed (prevents long range jumps)
% rf and cf - are the output row and column positions
% k-1 is the 'radius' of the final search matrix (due to pre-incrementing)

% to handle edge crashes and prevent premature stops at the edge edge of 
% the original nxn image we embedded it inside inside a (n+2k)x(n+2k)

% rs and cs - are the search row and column
rs = [];    % search row
k = 1;      % search radius
efail = 0;  % have we failed to find a point in the max search window?
while(isempty(rs))
    % make sure we haven't hit the edge
    % NOTE that this code fails if the image is not square
    if k<=maxk
        % pull out a test region around the last point
        subtest = bwimage(rp-k:rp+k,cp-k:cp+k);
        % find any non-zero points in the test region
        [rs,cs,vs]= find(subtest);
                
        % just in case we find nothing, let's increment the search
        % region parameter
        k = k+1;
    else
        % we must have hit maxk so we mark the point with 0's and get out
        rs=0;
        cs=0;
        % acknowledge our failure
        efail = 1;
    end
end

% assign the output values to the first point found
rf = rp+rs(1)-k;
cf = cp+cs(1)-k;
% deal with a failure
if efail
    rf = 0;
    cf = 0;
end
% reset the (post-indexed) search radius before returning it
k = k-1;

return
end
    