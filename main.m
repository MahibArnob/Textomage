%Read the image
RGBI = imread('../Images/img3.jpg');

%Convert to Gray Scale Image
GI = rgb2gray(RGBI);


figure
imshow(GI)
hold on
title('Grey Scale Image')
hold off

%Find threshold using otsuthresh for binarization
[counts,x] = imhist(GI,16);
%stem(x,counts) %User to plot graph
T = otsuthresh(counts);

%T1 = adaptthresh(GI, 0.9);
T1 = graythresh(GI);

%Convery gray scale image to binary image using the threshold
BI = imbinarize(GI,T);

figure
imshow(BI);
hold on
title('Binary Image')
hold off

%Find connected component labels
[LI] = bwlabel(BI);

%{
%Display connected components
for i = 1: num
    figure
    imshow((LI==i))
end
%}

stats = regionprops(LI, 'Area', 'BoundingBox', 'Eccentricity', ...
    'Solidity', 'Euler', 'FilledArea', 'Orientation');

%{
centroids = cat(1, stats.Centroid);

imshow(BI)
hold on
plot(centroids(:,1),centroids(:,2), 'b*')
hold off
%}

% Compute the aspect ratio using bounding box data.
bbox = vertcat(stats.BoundingBox);
%It return matrix with x, y, width, height.

w = bbox(:,3);
h = bbox(:,4);
aspectRatio = h./w;

% Threshold the data to determine which regions to remove. These thresholds
% may need to be tuned for other images.
idx = aspectRatio' > 0.96;
idx = idx & [stats.Eccentricity] > 0.32 ;
idx = idx & [stats.Solidity] < 0.66;
idx = idx & [stats.EulerNumber] < 1.000001 & [stats.EulerNumber] > -1.000001;
idx = idx & [stats.Area] > 101;
idx = idx & [stats.FilledArea] > 112;
idx = idx & [stats.Orientation] < 90 ;

keeperIdx = find(idx); %Returns indices of non-zero elements.

BI2 = ismember(LI, keeperIdx);
figure
imshow(BI2);

[LI2] = bwlabel(BI2);

stats = regionprops(LI2, 'BoundingBox');

% Compute the aspect ratio using bounding box data.
bbox = vertcat(stats.BoundingBox);

% Convert from the [x y width height] bounding box format to the [xmin ymin
% xmax ymax] format for convenience.
xmin = bbox(:,1);
ymin = bbox(:,2);
xmax = xmin + bbox(:,3) - 1;
ymax = ymin + bbox(:,4) - 1;

% Expand the bounding boxes by a small amount.
expansionAmount = 0.02;
xmin = (1-expansionAmount) * xmin;
ymin = (1-expansionAmount) * ymin;
xmax = (1+expansionAmount) * xmax;
ymax = (1+expansionAmount) * ymax;

% Clip the bounding boxes to be within the image bounds
xmin = max(xmin, 1);
ymin = max(ymin, 1);
xmax = min(xmax, size(LI2,2));
ymax = min(ymax, size(LI2,1));

% Show the expanded bounding boxes
expandedBBoxes = [xmin ymin xmax-xmin+1 ymax-ymin+1];
IExpandedBBoxes = insertShape(LI2,'Rectangle',expandedBBoxes,'LineWidth',2);

figure
imshow(IExpandedBBoxes)
title('Bounding Boxes Text')


% Compute the overlap ratio
overlapRatio = bboxOverlapRatio(expandedBBoxes, expandedBBoxes);

% Set the overlap ratio between a bounding box and itself to zero to
% simplify the graph representation.
n = size(overlapRatio,1);
overlapRatio(1:n+1:n^2) = 0;

% Create the graph
g = graph(overlapRatio);

% Find the connected text regions within the graph
componentIndices = conncomp(g);

% Merge the boxes based on the minimum and maximum dimensions.
xmin = accumarray(componentIndices', xmin, [], @min);
ymin = accumarray(componentIndices', ymin, [], @min);
xmax = accumarray(componentIndices', xmax, [], @max);
ymax = accumarray(componentIndices', ymax, [], @max);

% Compose the merged bounding boxes using the [x y width height] format.
textBBoxes = [xmin ymin xmax-xmin+1 ymax-ymin+1];

% Remove bounding boxes that only contain one text region
numRegionsInGroup = histcounts(componentIndices);
textBBoxes(numRegionsInGroup == 1, :) = [];

% Show the final text detection result.
ITextRegion = insertShape(LI2, 'Rectangle', textBBoxes,'LineWidth',2);

figure
imshow(ITextRegion)
title('Detected Text')

ocrtxt = ocr(LI2, textBBoxes);
[ocrtxt.Text]

%contour merging maybe