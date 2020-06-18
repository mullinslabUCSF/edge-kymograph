function kymo = generic_edge2(imgdir1,maskdir,resdir,name1,fname1,thresh,blur)
% generic_edge() - This function identifies the edge of a cell expressing a
% fluorescent protein. This edge becomes the ROI for identifying particles,
% such as VASP clusters. 
% RDM and KC 08/15/2016 
% revised RDM 9/12/2017
% revisited by RDM 5/2/2019
%
% imgdir - is the directory that contains the timelapse image sequence
% maskdir - is a directory that will hold edge-masked images for processing
% resdir - is the name of a directory to hold the processed images   
% forename - is the name prefix for reconstructing file names
% name - is a vector of name suffixes for reconstrucitng file names
% resultsdir - is a directory for saving the output (qCentroids)
% filename_tag - is a string used to identify the output files
% thresh - is the threshold used to discriminate objects from noise
% blur - is the degree of gaussian blurring to perform on initial image

% set the return directory to the current one
ret_dir = pwd;

% count the number of frames
n_fr = length(name1);
% initialize a cell array to hold the kymograph data
% the first position is the number of points (n)
% the second position is the nx2 matrix of x,y coordinates of the cell edge
% the third positoin is the nx1 vector of fluorescence intensities at edge
kymo = cell(n_fr,3);
dummy = zeros(n_fr,1);

% lwid - is the width of the region of interest at the cell edge
lwid = 20;
lrad = floor(lwid/2);
% perims - is a cell array holding all the 1-pixel-wide bw perimeters
perims = cell(lwid,1);
    
% grab the min and max intensities from the first image
cd(imgdir1)
I = imread([fname1,name1{1}]);
minI=min(min(I));
maxI=max(max(I));
maxI=maxI-(maxI-minI)/2;
% minI and maxI are the normalization parameters to use for the whole movie

% open a figure
%figure;

% run through all the frames in the directory and build the ROI's
for i=1:n_fr

    
    % go to the channel 1 image directory 
    cd(imgdir1)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    I = imread([fname1,name1{i}]);
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
    % troubleshooting
    %imshow(bw);
    %%%%%%%%%%%%%%%%%
    bw2 = imfill(bw,'holes');
    
    % may not need to use bwareaopen
    % but bw3 and bw4_perim are the two most useful bw images for analysis
    bw3 = bwareaopen(bw2,100); 
    bw4_perim = bwperim(bw3);
    
    hold on
    % dilates the detected perimeter to capture all the leading edge signal
    SE = strel('disk', lrad , 4);
    bw5 = imdilate(bw4_perim,SE);
    
    % use the edge mask (bw5) to mask the original image (I)
    I = I.*uint16(bw5);   
    
    % figure out where to break the perimeter for the kymograph
    if i==1
        % on the first pass let the user choose the boundary points
        [br_pt,end_pt,mid_pt]=boundary_set(bw4_perim);
        % now create figure 1 for overlaying cell boundaries from frames
        figure(1);
        % First we invert the image pixel values so that the background is white
        % and the cell perimeter is black. This gives better contrast for clicking
        % because the cursor is also black.
        invperim=imcomplement(bw4_perim);
        imshow(invperim)
    else
        % on subsequent passes choose boundary points close to the previous
        if ~bw4_perim(br_pt)
            [br_pt(1),br_pt(2),k]=findthread(bw4_perim,br_pt(1),br_pt(2),100);
        end
        if ~bw4_perim(end_pt)
            [end_pt(1),end_pt(2),k]=findthread(bw4_perim,end_pt(1),end_pt(2),100);
        end
        if ~bw4_perim(mid_pt)        
            [mid_pt(1),mid_pt(2),k]=findthread(bw4_perim,mid_pt(1),mid_pt(2),100);
        end
    end
    % order the points around the cell perimeter
    [path1,rgb_outline]=walkabout(bw4_perim,br_pt);
%     % add the new cell outline to the existing figure 1
%     figure(1);
%     hold on
%     plot(path1(:,2),path1(:,1),'color',[0 0 0]);
%     hold off
    
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
    path=path1(1:rend,:);
    num=length(path);
    % add the new ROI outline to the existing figure 1
    figure(1);
    cmap=parula(64);
    oin=floor((i/n_fr)*64);
    hold on
    plot(path1(:,2),path1(:,1),'color',[0 0 0]);
    plot(path(:,2),path(:,1),'color',cmap(oin,:));
    hold off
    
    % compute normal vectors at every point of the perimeter
    % normals = cell(num,1);
    [normals,RGBsl] = findnorm(bw4_perim,path,num,5,lrad);
    
    % load information into the kymograph cell array
    dummy(i) = num;
    kymo{i,1} = num;
    kymo{i,2} = path;
    % compute and store the intensity values from around the periphery
    kymo{i,3} = curvescan(I,normals,num);
    
    % use cross-correlation to align the intensity traces in the kymographs
    % (works, but ugly)
    if i>1
        % calculate the cross correlation
%         [r,lag] = xcorr(kymo{i-1,3},kymo{i,3});
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%         % subtract backgrounds
%         y1 = kymo{i-1,3} - mean(kymo{i-1,3});
%         y2 = kymo{i,3} - mean(kymo{i,3});
%         y1(y1<0) = 0;
%         y2(y2<0) = 0;
%         [r,lag] = xcorr(y1,y2);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % smooth out the cross correlation a bit
%         windowSize = 31;
%         fb = (1/windowSize)*ones(1,windowSize);
%         fa = 1;
%         r = filter(fb,fa,r);
%           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % identify the offset by finding the max of the cross correlation
%         tl = lag(find(r==max(r)));
%         % shift the data and path
%         % old way
%         %kymo{i,2} = circshift(kymo{i,2},tl);
%         %kymo{i,3} = circshift(kymo{i,3},tl);
%         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % now try a second round of cross correlation without bkgd subtract
        % calculate the cross correlation
        [r,lag] = xcorr(kymo{i-1,3},kymo{i,3});
        tl = lag(find(r==max(r)));
        % shift the data and path
        kymo{i,2} = circshift(kymo{i,2},tl);
        kymo{i,3} = circshift(kymo{i,3},tl);
        % limit the damage of the circular perutation
        if tl>0
            kymo{i,3}(1:tl)=0;
        end
    end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % display the computed cell outline
%         subplot(1,3,2);
        figure(2);
        imshow(RGBsl);
        % mark the boundary points
        hold on
        rectangle('Position',[br_pt(2)-10 br_pt(1)-10 20 20],'Curvature',[1 1],'EdgeColor','cyan');
        rectangle('Position',[mid_pt(2)-10 mid_pt(1)-10 20 20],'Curvature',[1 1],'EdgeColor','cyan');
        rectangle('Position',[end_pt(2)-10 end_pt(1)-10 20 20],'Curvature',[1 1],'EdgeColor','cyan');
        plot(path(:,2),path(:,1),'ow');
        hold off
        %imshow(rgb_outline);
        title('cell outline');
%         subplot(1,3,3);
        % display a sample slope calculation
%         subplot(1,3,3);
        % display the bw5-masked image        
%         imshow(I,[minI maxI]); 
%         title(['masked with bw5']);

        % now we save the images for making a movie
        % first, grab the image from the figure
        F = getframe(gcf);
        % go to the results directory
        cd(resdir)
    
        % save the file with a good filename
        filename = ['edge_' name1{i} '.jpg'];
        imwrite(F.cdata,filename);
        hold off
        
    
    % now save the 16-bit masked image for further processing
    % go to the dedicated 'masked image directory'
    cd(maskdir)
    
    % make a new file name
    filename = ['masked_' name1{i}];
    imwrite(I,filename);
end

% % add a color bar to figure(1)
% figure(1)
% hold on
% colormap(1,parula);
% colorbar('south')
% hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build a kymograph image.
stretch = 5;
width = max(dummy);
im_kym = zeros(width,stretch*n_fr);
% mark the corner of the kymograph

for i=1:n_fr
    % the commented code centers the lines of the kymo which is a mistake
    %init = 1+floor((width-kymo{i,1})/2);
    for j=1:stretch
        im_kym(1:kymo{i,1},1+(i-1)*stretch+j-1) = kymo{i,3};
    end
end
hold off    % just in case
figure(3);
cmap = colormap(hot(128));
imshow(im_kym,[0 max(max(im_kym))],'colormap',jet(128));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% there's no place like home
cd(ret_dir)

return
end

function [path,rgb_outline]=walkabout(bwimage,start)
% this function takes a binary image of a closed (Jordan) curve and
% returns an ordered array of x,y indices that walks, stepwise around the
% entire curve

% path - is a two-column matrix containing the ordered x and y coordinates
% num - is the number of points in the closed curve
% rgb_outline - is a color RGB version of the cell perimeter color coded 
% bwimage - is the input image containing the one-pixel-wide curve
% start - is a pair of coordinates defining the starting position
%           N.B. if start = [0 0] the alogrithm automatcally breaks
%           symmetry at a point near the middle left edge of the cell

% first find all of the non-zero elements in the image (row and column are
% the indices of the non-zero elements; v is the value of each pixel)
[row,col,v]=find(bwimage);

% compute the length of the output matrix
num = length(row);
% here is a backstop variable to hold the output matrix length
% this code needs to be cleaned up because ultimate arc length is ambig.
early_stop = 0;

% initialize vectors to hold the output (sorted) rows and columns
rf = zeros(num,1);
cf = zeros(num,1);

if start(1)==0
    % here we break symmetry by picking the mid-left edge of the cell
    [r1,v]=find(row==260);
    cp=min(col(r1));
    rp=260;
else
    % here, we use the input starting point
    rp = start(1);
    cp = start(2);
end

% load up the first position
rf(1) = rp;
cf(1) = cp;

% now run through the connected points on the curve and order them
for i=1:num-1
    % remove the point at the current position
    if i==1
        bwimage(rp-12:rp+12,cp-12:cp+12)=zeros(25);
    else
        bwimage(rp,cp) = 0;  
    end
    % use findthread() to get the next point
    [rp,cp,k]=findthread(bwimage,rp,cp,20);
    rf(i+1)=rp;
    cf(i+1)=cp;
end

% % make sure that we walked counterclockwise around the edge of the cell
% if sum(rf(1:floor(num/2)))<sum(rf(ceil(num/2):num))
%     rf = flip(rf);
%     cf = flip(cf);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the output parameters
% take account of the fact that the trace probably ends early
if early_stop == 0
    early_stop = length(rf);
end
path = [rf(1:early_stop) cf(1:early_stop)];
num = early_stop;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create a color-coded image of the cell outline
cmap = hot(128);
rgb_outline=zeros(512,512,3);
for i=1:num
    rgb_outline(path(i,1),path(i,2),:)=cmap(ceil(i*(98/num)+30),:);
end

return
end

function [norm,RGB_slope] = findnorm(perimIM,path,npts,smrad,outrad)
% this function uses the matrix 'path' to navigate around the perimeter
% image 'perimIM' and calculate normal vectors to every point on the
% perimeter. These are returned in the cell structure 'norm'
%
% perimIM
% path
% npts
% smrad
% outrad

% position 1 holds the columns (y-data) for the kymograph evaluation
% position 2 holds the rows (x-data) for the kymograph evaluation
norm = cell(npts,3);

% run around the perimeter
for i=1:npts
    % set the current point
    xo = path(i,2);
    yo = path(i,1);
    
    % pull out a submatrix centered on (xo, yo)
    subtest = perimIM(yo-smrad:yo+smrad,xo-smrad:xo+smrad);
    % identify the perimeter points within the submatrix
    [ys,xs,vs]= find(subtest);
    % center the submatrix points on the origin
    xs = xs-(smrad+1);
    ys = ys-(smrad+1);
    
    % use linear regression to find the best-fit slope of a line that runs
    % through the origin: i.e. y = mx OR x = my
    % to deal with large slopes, we will calculate slope both ways
    
    % since y = mx, the 'estimator matrix' is just xs
    % the observation vector is ys
    % the results are contained in res
    
    % solve for the best fit coefficient with x as the independent variable
    xres = (xs'*xs)\xs'*ys;
    % solve for best fit coefficient with y as the independent variable
    yres = (ys'*ys)\ys'*xs;
    
    % N.B. To find a normal to the tangent we just invert the slope and
    % multiply by -1. The trick here is that the y-based slope is the
    % inverse of the x-based slope
    
    % choose the smaller of the two slopes to calculate the points that lie
    % on the normal line: this avoids problems with infinite slopes:
    % if the x-based slope is larger, then use the y-based slope...
    if abs(xres) > abs(yres) || isnan(xres)
        % do the calculation right
        yso = [-outrad:outrad]';
        xso = round(yres*yso);
        % then rotate the result by 90 degrees to create a normal
        norm{i,1} = -xso + yo;
        norm{i,2} = yso + xo;
    % similarly, if the y-based slope is larger...    
    else
        % do the calculation right
        xso = [-outrad:outrad]';
        yso = round(xres*xso);
        % then rotate the result by 90 degrees to create a normal
        norm{i,1} = -xso + yo;
        norm{i,2} = yso + xo;
    end
    % by using the negative of the reciprocal slopes we automatically
    % rotate the tangent vector by 90 degrees to create the normal vector
    
    % finally, save the number of points in the normal vector
    norm{i,3} = length(norm{i,1});
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct an RGB image containing the color-coded slopes
% initialize the RGB image
RGB_slope = zeros(512,512,3);
cmap = hot(128);
for i=1:npts 
    % load the output RGB image
    %if i<180 || i>220    % this is to add a polarity mark for debugging
    ind = 1+floor(127*abs(i/npts));
    for j=1:norm{i,3}
%         norm{i,1}(j)
%         norm{i,2}(j)
        RGB_slope(norm{i,1}(j),norm{i,2}(j),:) = cmap(ind,:);
    end
    %end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return
end

function intens = curvescan(image,norms,npts)
% this function computes the maximum intensity values along a set of normal
% vectors defined on a grayscale image.

% intens - is the output vector containing the measured intensities
% image - is the original 2D imaging data
% norms - is a cell array containing x and y values of normals at each
%         point on the curve 
% npts - is the number of points in the perimeter line

% initialize the output parameter
intens = zeros(npts,1);

% run around the curve and compute the intensities
for i=1:npts
    vals = zeros(norms{i,3},1);
    for j=1:norms{i,3}
        vals(j) = image(norms{i,1}(j),norms{i,2}(j));
    end
    intens(i) = max(vals);
end

return
end

function [startpt,endpt,midpt]=boundary_set(bwimage)
% Function boundary_set enables the user to define a region of interest 
% along the perimeter of the input image(bwimage) using the mouse. 
% The function returns the coordinates of the start and end points together
% with a midpoint required to define which segment of the perimeter to
% analyze.

figure;
% First we invert the image pixel values so that the background is white
% and the cell perimeter is black. This gives better contrast for clicking
% because the cursor is also black.
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

    % check to see whether we have hit the perimeter image on the first try
    if ~bwimage(y(i),x(i))
        % find the point in the perimeter image closest to the chosen pixel
        % to do this we will use the function findthread()
        [y(i),x(i),dumk]=findthread(bwimage,y(i),x(i),100);  
    end
    
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

% rp and cp - are the input row and coumn position for the search
% bwimage - is a logical image containing the perimeter of an object
% maxk - is the maximum search radius allowed (prevents long range jumps)
% rf and cf - are the output row and column positions
% k-1 is the 'radius' of the final search matrix 

% rs and cs - are the search row and column
rs = [];    % search row
k = 1;      % search radius
while(isempty(rs))
    % make sure we haven't hit the edge
    if rp-k>0 && cp-k>0 && rp+k<length(bwimage) && cp+k<length(bwimage) && k<=maxk
        % pull out a test region around the last point
        subtest = bwimage(rp-k:rp+k,cp-k:cp+k);
        % find any non-zero points in the test region
        [rs,cs,vs]= find(subtest);
                
        % just in case we find nothing, let's increment the search
        % region parameter
        k = k+1;
    else
        rs=0;
        cs=0;
        k=0;
    end
end

% assign the output values to the first point found
rf = rp+rs(1)-k;
cf = cp+cs(1)-k;
% reset the (post-indexed) search radius before returning it
k = k-1;

return
end
