function varargout = GUI(varargin)
% GUI MATLAB code for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 23-Mar-2017 16:15:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Read the image

[filename,pathname]=uigetfile({'*.*'},'Select file');
var=strcat(pathname,filename);
RGBI=imread(var);

figure
imshow(RGBI)
hold on
title('Original Image')
hold off

%Convert to Gray Scale Image
GI = rgb2gray(RGBI);

%Find threshold using otsuthresh for binarization
[counts] = imhist(GI,16);
T = otsuthresh(counts);

%Convery gray scale image to binary image using the threshold
BI = imbinarize(GI,T);

ones = numel(find(BI));
total = numel(BI);

if ones > (total-ones)
    BI = ~BI;
end

%Find connected component labels
[LI] = bwlabel(BI);

stats = regionprops(LI, 'Area', 'BoundingBox', 'Eccentricity', ...
    'Solidity', 'Euler', 'FilledArea', 'Orientation', 'Perimeter');

% Compute the aspect ratio using bounding box data.
bbox = vertcat(stats.BoundingBox);
%It return matrix with x, y, width, height.

w = bbox(:,3);
h = bbox(:,4);
aspectRatio = h./w;

% Threshold the data to determine which regions to remove. These thresholds
% may need to be tuned for other images.
idx = aspectRatio' > 0.74; %0.96
idx = idx & [stats.Eccentricity] > 0.20 ; %0.32
idx = idx & [stats.Solidity] < 1.00001; %0.66
idx = idx & [stats.EulerNumber] < 1.000001 & [stats.EulerNumber] > -1.000001;
idx = idx & [stats.Area] > 101;
idx = idx & [stats.FilledArea] > 112;
idx = idx & [stats.Orientation] < 90.0001 ; %90
idx = idx & ([stats.Perimeter] <= 80.7 | [stats.Perimeter] >= 80.8);
idx = idx & ([stats.Perimeter] <= 143.15 | [stats.Perimeter] >= 143.17);

keeperIdx = find(idx); %Returns indices of non-zero elements.

BI2 = ismember(LI, keeperIdx);

ones = numel(find(BI2));
total = numel(BI2);

if ones > (total-ones)
    BI2 = ~BI2;
end

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
expansionAmount = 0.1;
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

%Merge the columns to make a matrix.
expandedBBoxes = [xmin ymin xmax-xmin+1 ymax-ymin+1];

BBI = insertShape(LI2,'Rectangle',expandedBBoxes,'LineWidth',2);

% Compute the overlap ratio
overlapRatio = bboxOverlapRatio(expandedBBoxes, expandedBBoxes);

% Set the overlap ratio between a bounding box and itself to zero to
% simplify the graph representation.
n = size(overlapRatio,1);
overlapRatio(1:n+1:n^2) = 0;

% Create the graph
g = graph(overlapRatio);
plot(g)

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

ocrtxt = ocr(LI2, textBBoxes);
[ocrtxt.Text]

set(handles.text2, 'String', ocrtxt.Text);
