function [startline,endline,midline,bkgd]=boundary_lines(bwimage)
% Function boundary_lines enables the user to define a region of interest 
% along the perimeter of the input image(bwimage) using the mouse. 
% The function returns the coordinates of the start and end points together
% with a midpoint required to define which segment of the perimeter to
% analyze.

% startline - line defining the starting position of ROIs 
% endline - line defining the ending position of leading edge ROIs
% midline - line defining (roughly) the midpoint of the ROIs
% bkgd - rectangle defining region used for background subtraction

% First we display the input image with a few instructions
figure(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bwimage2=imcomplement(bwimage);
bwimage2=bwimage;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
text_str = 'Click two points each to define ROI: (1)start, (2)end, (3)midpoint';
imshow(bwimage2)
text(15,15,text_str);

% initialize the data values
x = zeros(8,1);
y = zeros(8,1);

% set up a colormap for plotting the breakpoints
cmap=hot(8);
hold on
% now get six points from the image
% three pairs of images
ind = 0;
for i=1:3
    for j=1:2
        % pre increment
        ind = ind+1;
        [x(ind),y(ind)]=ginput(1);
        % make sure the values are well behaved
        y(ind)=floor(y(ind));
        x(ind)=floor(x(ind));

%     % check to see whether we have hit the perimeter image on the first try
%     if ~bwimage(y(i),x(i))
%         % find the point in the perimeter image closest to the chosen pixel
%         % to do this we will use the function findthread()
%         [y(i),x(i),dumk]=findthread(bwimage,y(i),x(i),100);  
%     end
%     
        % now plot the chosen point on the image
        cind = ceil(ind/2);       % color index
        rectangle('Position',[x(ind)-5 y(ind)-5 10 10],'Curvature',[1 1],'EdgeColor',cmap(cind,:));
    end
    % now connect the dots to make a line
    line([x(ind-1) x(ind)],[y(ind-1) y(ind)],'color','cyan');
end

% now identify a rectangle used for background subtraction
for i=1:2
    [x(6+i),y(6+i)]=ginput(1);
    % make sure the values are well behaved
    y(6+i)=floor(y(6+i));
    x(6+i)=floor(x(6+i));
    % plot the corner points of the background rectangle
    rectangle('Position',[x(6+i)-5 y(6+i)-5 10 10],'Curvature',[1 1],'EdgeColor',cmap(cind,:));
end

% connect the corners to make a rectangle
line([x(7) x(8)],[y(7) y(7)],'color','cyan');
line([x(7) x(8)],[y(8) y(8)],'color','cyan');
line([x(7) x(7)],[y(7) y(8)],'color','cyan');
line([x(8) x(8)],[y(7) y(8)],'color','cyan');

hold off

% load up the return variables
startline=[y(1),x(1);y(2),x(2)];
endline=[y(3),x(3);y(4),x(4)];
midline=[y(5),x(5);y(6),x(6)];
bkgd=[y(7),x(7);y(8),x(8)];
% because in matlab, [col row] = [y,x]

return
end

