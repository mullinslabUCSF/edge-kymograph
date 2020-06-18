function outimage = makeline(xy1,xy2,inpimage)
% line_img - creates a logical image 'line_img' with a line segment that
% runs from xy1 to xy2
% 
% xy1 and xy2 - are the endpoints for the line
% inpimage is a reference image with the proper dimensions for the output
% line_img - is an output bw image with the same dimensions as inpimage
% with the line segment drawn in

deltax = xy2(1)-xy1(1);
deltay = xy2(2)-xy1(2);

slope=abs(double(deltay)/double(deltax));
invsl=abs(double(deltax)/double(deltay));

% make a submatrix for computation
submx = zeros(abs(deltax)+1,abs(deltay)+1);
% light up the origin because it is, by definition, on the line segment
submx(1,1)=1;
% run through the x values and light up corresponding y's
for i = 1:abs(deltax)
    yin=round(slope*double(i)+1);
    submx(i+1,yin)=1;
end
% run through the y values and light up corresponding x's
for i = 1:abs(deltay)
    xin=round(invsl*double(i)+1);
    submx(xin,i+1)=1;
end

% now that we have drawn the line, let's orient it correctly
% if xy1(1)<xy2(1) AND xy1(2)<xy2(2) do nothing!
% or if xy1(1)>xy2(1) AND xy1(2)>xy2(2) do nothing!
if xy1(1)<xy2(1) && xy1(2)>xy2(2)
    submx=flipud(submx);
elseif xy1(1)>xy2(1) && xy1(2)<xy2(2)
    submx=fliplr(submx);
end

% next, load the line segment into the input image to create the output
% image
x1 = min(xy1(1),xy2(1));
x2 = max(xy1(1),xy2(1));
y1 = min(xy1(2),xy2(2));
y2 = max(xy1(2),xy2(2));
outimage=inpimage.*0;
outimage(x1:x2,y1:y2) = submx;
outimage=logical(outimage);

% finally, expand the line a bit to make sure it scores intersection hits
SE = strel('disk', 2, 4);
outimage = imdilate(outimage,SE);

return
end

