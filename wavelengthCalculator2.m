%Load in several vertical slices of data, call
%Pinch_Location_Calculator(), then curve fit all the X or Y locations to
%a polynomial curve to determine wavelength. 

%Clear/close all, manually input axial locations and frame time/shot to analyze
clear all
clc
close all

%Currently working on throw/catch commands to handle incorrect location
%inputs

locations = [367,387,407,427,447,487,507,527,547,567]; %Location of the pixels grabbed from Kirana (X=437,447,...)
%locations = [437,447,457,477,487,497,507,517,527,547];
%locations = [447,487,507,527];
%locations = [275,295,315,335,375,395,415,435];
%locations = [315,375,395,415,435,455];
%locations = [325,345,385,405,425,445,465,485];
shotInput = 210825008;
timeInput = 50.75E-6; %sec

locInput = 'P15';
locInputPixel = 467;
%locInputPixel = 365;

windowInput1 = 88; %Top window limits for Aug 25
windowInput2 = 212;
windowInput3 = 535; %Bottom window limits
windowInput4 = 674;

% windowInput1 = 32; %Top window limits for Sept 22
% windowInput2 = 165;
% windowInput3 = 511; %Bottom window limits
% windowInput4 = 663;

% windowInput1 = 38; %Top window limits for Oct 14
% windowInput2 = 170;
% windowInput3 = 510; %Bottom window limits
% windowInput4 = 660;

%Analyze both windows at once, so 'windowSelection' is not needed
%windowSelection = 2; %2 for bottom window (Y), 1 for top window (X)

%Analyze bottom plot -- need Y location of visual
%For loop reads in all values near P15

counter = 1;
locationsErrorCounter = 1;
for j = 1:length(locations)
%     if j>length(locations)
%         break
%     end
    [visual, ~, ~, ~,peakCenterMatrix] = Pinch_Location_Calculator(shotInput,timeInput,num2str(locations(j)),windowInput1,windowInput2,windowInput3,windowInput4,false);
    %These two matrices hold the X/Y pinch locations from visual data
    %analysis
    %X & Y Loc are the locations in centimeters from the center
    XLoc(j) = visual(1); %2 for bottom window (Y), 1 for top window (X)
    YLoc(j) = visual(2);
    
    %X & YLocPixels are the locations in pixels on the picture, these
    %variables are used for plotting points over Kirana picture.
    XLocPixels(j) = peakCenterMatrix(1);
    YLocPixels(j) = peakCenterMatrix(2);
    
%     if sum(peakCenterMatrix)~=0
%        XLoc(counter) = visual(1); %2 for bottom window (Y), 1 for top window (X)
%        YLoc(counter) = visual(2);
%        XLocPixels(counter) = peakCenterMatrix(1);
%        YLocPixels(counter) = peakCenterMatrix(2); 
%        counter = counter+1;
%     else
%         locationsError(locationsErrorCounter) = j;
%         locationsErrorCounter = locationsErrorCounter + 1;
%         %locations(locations==locations(j))=[];
%     end
end
%Repeat pinch location calculation w/ P15
%Set last argument to true in order to compare the visual data to jScope
%data
[visual, t, tplus1, tminus1,peakCenterMatrix] = Pinch_Location_Calculator(shotInput,timeInput,locInput,windowInput1,windowInput2,windowInput3,windowInput4,true);
XLocP15 = visual(1);
YLocP15 = visual(2);
if length(peakCenterMatrix) == 1 %this if-else statement handles an cases where P15 does not contain pinches for both windows
    XLocPixelsP15 = 0;
    YLocPixelsP15 = peakCenterMatrix(1);
else
    XLocPixelsP15 = peakCenterMatrix(1);
    YLocPixelsP15 = peakCenterMatrix(2);
end

%This for loop will determine where to insert P15/other desired location
%into the rest of the pixel locations in the location/X/Y matrix
for j = 1:length(locations)
    if locations(j)<locInputPixel && locations(j+1)>locInputPixel
        locIndex = j;
    end
end
%disp(locIndex)

%This section of code inserts P15 into the matrices as described above. We
%are assuming that the desired location for analysis not the lowest or
%highest pixel in the matrix.
locations = [locations(1:locIndex),locInputPixel,locations(locIndex+1:end)];
XLoc = [XLoc(1:locIndex), XLocP15, XLoc(locIndex+1:end)];
YLoc = [YLoc(1:locIndex), YLocP15, YLoc(locIndex+1:end)];
XLocPixels = [XLocPixels(1:locIndex), XLocPixelsP15, XLocPixels(locIndex+1:end)]; %Added in these two last lines to try and plot everything in pixels only, and not in centimeters
YLocPixels = [YLocPixels(1:locIndex), YLocPixelsP15, YLocPixels(locIndex+1:end)];
% locations = sort([locations,locInputPixel]);
% XLoc = sort([XLoc,XLocP15]);
% YLoc = sort([YLoc,YLocP15]);
% XLocPixels = sort([XLocPixels,XLocPixelsP15]);
% YLocPixels = sort([YLocPixels,YLocPixelsP15]);

%% Sine wave curve fitting
%Old code encased in this for loop (curve fitting to sine wave)
for d = 1
% %% Sine wave curvefit in X direction
% y = XLoc;
% x = locations;
% 
% yu = max(y);
% yl = min(y);
% yr = (yu-yl);                               % Range of ‘y’
% yz = y-yu+(yr/2);
% zx = x(yz .* circshift(yz,[0 1]) <= 0);     % Find zero-crossings
% per = 2*mean(diff(zx));                     % Estimate period
% ym = mean(y);                               % Estimate offset
% fit = @(b,x)  b(1).*(sin(2*pi*x./b(2) + 2*pi/b(3))) + b(4);    % Function to fit
% fcn = @(b) sum((fit(b,x) - y).^2);                              % Least-Squares cost function
% options = optimset('MaxFunEvals',10000);
% s = fminsearch(fcn, [yr;  per;  -1;  ym],options)                       % Minimise Least-Squares
% xp = linspace(min(x),max(x));
% figure(6)
% plot(x,y,'b',  xp,fit(s,xp), 'r')
% title('Top Window, X Locations')
% grid
% 
% % Calculate wavelength in cm using s(2) in the X direction
% wavelengthCM_X = s(2)/(windowInput2-windowInput1)*5;
% disp(wavelengthCM_X)
% %% Sine wave curvefit in Y direction
% clear x y yu y1 yr yz zx per ym fit fcn s xp
% y = YLoc;
% x = locations;
% 
% yu = max(y);
% yl = min(y);
% yr = (yu-yl);                               % Range of ‘y’
% yz = y-yu+(yr/2);
% zx = x(yz .* circshift(yz,[0 1]) <= 0);     % Find zero-crossings
% per = 2*mean(diff(zx));                     % Estimate period
% ym = mean(y);                               % Estimate offset
% fit = @(b,x)  b(1).*(sin(2*pi*x./b(2) + 2*pi/b(3))) + b(4);    % Function to fit
% fcn = @(b) sum((fit(b,x) - y).^2);                              % Least-Squares cost function
% s = fminsearch(fcn, [yr;  per;  -1;  ym],options)                       % Minimise Least-Squares
% xp = linspace(min(x),max(x));
% figure(7)
% plot(x,y,'b',  xp,fit(s,xp), 'r')
% title('Bottom window, Y Locations')
% grid
% 
% % Calculate wavelength in cm using s(2) in the Y direction
% wavelengthCM_Y = s(2)/(windowInput4-windowInput3)*5;
% disp(wavelengthCM_Y)

% %% Calculate effective wavelength
% kx = 2*pi/wavelengthCM_X;
% ky = 2*pi/wavelengthCM_Y;
% 
% effWavelength = sqrt(kx^2+ky^2);

%This data should be received from Pinch_Location_Calculator(), no need to
%call data from MDSplus
%May change Pinch_Location_Calculator to call jScope instead of local files
%to be more robust
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate displacement of visual and jScope data and plot against
%effWavelength/kx/ky?
% mdsconnect('zappa.zap');
% mdsopen('zaphd',shotInput);
% MagTime = 1e6*mdsvalue('dim_of(\b_p0_0_t)');
% XLoc_probes = 100*mdsvalue('.5 * \R_OUTER * SIGMUL(\M_1_P15_NORM, COS(\PHI_M_1_P15 - $PI / 8))');
% YLoc_probes = 100*mdsvalue('.5 * \R_OUTER * SIGMUL(\M_1_P15_NORM, SIN(\PHI_M_1_P15 - $PI / 8))');
% mdsclose;
% mdsdisconnect;

% Currently have time/x/y locations stored in .mat files, need to convert
% to txt files and place in respective folders so that we can compare to
% visual data if MDSplus doesn't work
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%% Curve fit pinch to a polynomial and calculate radius of curvature

syms x y
%Curve fit the original points (pre-rotation), assume 3rd order polynomial
polyXpre = polyfit(locations, XLocPixels,3);
fxpre(x) = polyXpre(1)*x^3 + polyXpre(2)*x^2 + polyXpre(3)*x + polyXpre(4);
polyYpre = polyfit(locations, YLocPixels,3);
fypre(y) = polyYpre(1)*y^3 + polyYpre(2)*y^2 + polyYpre(3)*y + polyYpre(4);

%This block of code will attempt to rotate the plot to account for a linear
%bias
slopeX = (XLocPixels(1)-XLocPixels(end))/(locations(1)-locations(end));
angleX = -atan(vpa(slopeX));
slopeY = (YLocPixels(1)-YLocPixels(end))/(locations(1)-locations(end));
angleY = -atan(vpa(slopeY));
for j = 1:length(locations)
    %These two matrices contain the new rotated locations of the centroid
    %in the format of mat = [new location pixel, new X or Y location pixel;
    %...]
    
    %Equation we're using is [x';y'] = [cos(th) -sin(th); sin(th) cos(th)]*[x; y], where x is from the locations matrix and y is from
    %the X or YLocMatrix variable
    
    %Since we're rotating around P15, we need to plug in [x-x_P15;
    %y-y_P15], and then add [x_P15; y_P15] at the end
    xmat(j,1:2)=([cos(angleX) -sin(angleX); sin(angleX) cos(angleX)]*[(locations(j)-locations(locIndex+1));(XLocPixels(j)-XLocPixels(locIndex+1))]+[locations(locIndex+1);XLocPixels(locIndex+1)])';
    ymat(j,1:2)=([cos(angleY) -sin(angleY); sin(angleY) cos(angleY)]*[(locations(j)-locations(locIndex+1));(YLocPixels(j)-YLocPixels(locIndex+1))]+[locations(locIndex+1);YLocPixels(locIndex+1)])';
end
%Curve fit using new rotated pixels
polyX = polyfit(double(xmat(:,1)), double(xmat(:,2)),3);
fx(x) = polyX(1)*x^3 + polyX(2)*x^2 + polyX(3)*x + polyX(4);
polyY = polyfit(double(ymat(:,1)), double(ymat(:,2)),3);
fy(y) = polyY(1)*y^3 + polyY(2)*y^2 + polyY(3)*y + polyY(4);

%calculate magnitude of curvature vector, K
Kx(x) = abs(diff(fx(x),2))/(1+(diff(fx(x))^2))^(3/2);
Ky(y) = abs(diff(fy(y),2))/(1+(diff(fy(y))^2))^(3/2);
for j = 1:length(locations)
    KxMat(j) = Kx(locations(j));
    KyMat(j) = Ky(locations(j));
end

%Calculate avg curvature 
KxAvg = double(mean(KxMat));
KyAvg = double(mean(KyMat));

%Determine amplitude of curve by doing (max+min)/2
amplitudeX = double((max(xmat(:,2))-min(xmat(:,2)))/2);
amplitudeY = double((max(ymat(:,2))-min(ymat(:,2)))/2);

cX = sqrt(max(KxMat)/amplitudeX);
cY = sqrt(max(KyMat)/amplitudeY);

wavelengthX = 2*pi/cX;
wavelengthY = 2*pi/cY;

wavelengthX_CM = double(wavelengthX/(windowInput2-windowInput1)*5);
wavelengthY_CM = double(wavelengthY/(windowInput4-windowInput3)*5);
%Attempting this block of code but slightly edited to show both curves at
%once and with pixels as the y axis (instead of cm)
figure(6)
I1 = imread(strcat(num2str(shotInput),'/',num2str(shotInput),' Picture.jpg'));
minx = 0;
maxx = 924;
miny = 0;
maxy = 768;
imagesc([minx maxx],[miny maxy],I1)
hold on
scatter(locations,XLocPixels,50,'filled','r')
scatter(locations,YLocPixels,50,'filled','r')
plot(locations(1):locations(end),fxpre(locations(1):locations(end)),'k','LineWidth',2)%Plotting original curve fits
plot(locations(1):locations(end),fypre(locations(1):locations(end)),'k','LineWidth',2)
%plot(xmat(1,1):xmat(end,1),fx(xmat(1,1):xmat(end,1)),'LineWidth',2)%Plotting rotated curve fits
%plot(ymat(1,1):ymat(end,1),fy(ymat(1,1):ymat(end,1)),'LineWidth',2)
%text(400,200,strcat('Wavelength: ',num2str(wavelengthX_CM),'cm'),'FontSize',14,'Color','r')
text(250,200,strcat('X Axis Average Curvature: ',num2str(KxAvg)),'FontSize',14,'Color','r')
text(250,250,strcat('Wavelength: ',num2str(wavelengthX_CM),'cm'),'FontSize',14,'Color','r')
text(250,450,strcat('Wavelength: ',num2str(wavelengthY_CM),'cm'),'FontSize',14,'Color','r')
text(250,400,strcat('Y Axis Average Curvature: ',num2str(KyAvg)),'FontSize',14,'Color','r')

hold off

%Old plotting code below, not currently used
for d=1
% figure(6)
% I1 = imread(strcat(num2str(shotInput),'/',num2str(shotInput),' Picture Top.jpg'));
% %I1 = imread('210825008/210825008 Picture 2 Edited Top.jpg');
% minx = 238;
% maxx = 612;
% miny = 2.5;
% maxy = -2.5;
% imagesc([minx maxx],[miny maxy],I1)
% set(gca,'YDir','normal')
% hold on
% scatter(locations,XLoc,50,'filled')
% plot(locations(1):locations(end),fx(locations(1):locations(end)))
% %axis([335,693,-2.5,2.5])
% hold off
% 
% figure(7)
% I2 = imread(strcat(num2str(shotInput),'/',num2str(shotInput),' Picture Bottom.jpg'));
% %I2 = imread('210825008/210825008 Picture 2 Edited Bottom.jpg');
% minx = 107;
% maxx = 657;
% miny = 2.5;
% maxy = -2.5;
% imagesc([minx maxx],[miny maxy],I2)
% set(gca,'YDir','normal')
% hold on
% scatter(locations,YLoc,50,'filled')
% plot(locations(1):locations(end),fy(locations(1):locations(end)))
% %axis([8,726,-2.5,2.5])
% hold off
end