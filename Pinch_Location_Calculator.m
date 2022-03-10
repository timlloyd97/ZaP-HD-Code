%Purpose: function that will accept .csv files from Kirana and will compare them to
%.jpg files of frames from Kirana. 

%INPUTS:
%shotInput is in YYMMDD### for shot number ###

%timeInput should be in seconds i.e 15E-6 sec

%locInput should be a pixel number of where you are analyzing the Kirana
%frame i.e 400. Currently code works when .csv files are analyzing a column
%of pixels

%windowInput1 is the top y pixel location of the top-down view window
%windowInput2 is the bottom y pixel location of the top-down view of the
%window
%windowInput3 is the top y pixel location of the side view of the window
%windowInput1 is the bottom y pixel location of the side view of the
%window

%compareToProbe is either true/false. Setting to true will attempt to pull
%centroid location from jScope to compare against Kirana visual analysis.
%Setting to false will simply analyze .csv and .jpg files given to it

%Naming convention: currently, each shot should have its own folder with
%the naming convention YYMMDD###/YYMMDD### Pixel# SS_SS.csv, where
%Pixel# is what pixel column you're analyzing, and SS_SS is the time of the
%frame in microseconds i.e 50.75us = 55_75. So shot 210825008 for pixel
%column 407 at 50.75us should have naming convention 210825008/210825008
%407 50_75.csv

%OUTPUTS:
%visual: gives the x and y location of the centroid, according to the
%Kirana frame in cm, relative to the z-axis. Given in array form
%visual=[x,y]

%t: gives x and y location of centroid at the same time you're analyzing,
%according to the B Field probes from jScope. Will return an actual value
%if compareToProbes=True, otherwie returns [0,0]

%tplus1: gives x and y location of centroid at exactly 1us later than the time you're analyzing,
%according to the B Field probes from jScope. Will return an actual value
%if compareToProbes=True, otherwie returns [0,0]

%tminus1: gives x and y location of centroid at exactly 1us earlier than the time you're analyzing,
%according to the B Field probes from jScope. Will return an actual value
%if compareToProbes=True, otherwie returns [0,0]

%peakCenterMatrix: returns pixel locations of where the pinch is located.
function [visual, t, tplus1, tminus1,peakCenterMatrix] = Pinch_Location_Calculator(shotInput,timeInput,locInput,windowInput1,windowInput2,windowInput3,windowInput4,compareToProbe)
%% Calculate pinch location
% For Aug 04, X=562 is P15 (at least for last few shots, disregard this
% note for previous shots) %%%%%%%%%%%%%%%%%%%
% For Aug 11, X=526 is P15 for all shots

% clear all
close all
% clc

%A is a column matrix with rows of the y pixels from the kirana camera.
%A higher y indicates a lower pixel on the Kirana viewing screen 
%[(X,Y) = (0,0) is the top left]

%anylyzing top, bottom, neither?
%wavelengthPos = 1 for top, 2 for bottom, 0 for neither
% wavelengthPos = 0;
% shot = 210811015;
% time = 43.68e-6; %sec
%wavelengthPos = wavePos;
shot = shotInput;
time = round(timeInput,10);

timeSec = floor(time*1e6);
timeDecSec = erase(num2str(time*1e6-timeSec),'0.');
%timeDecSec = 100*(time*1e6-timeSec);

%Read in visual data from Kirana csv file
if length(num2str(timeDecSec)) == 1 %This if statement prevents errors in case time is a whole number (Ex: 20us)
   fileName = [num2str(shot),'\',num2str(shot), ' ', locInput,' ',num2str(timeSec),'_',num2str(timeDecSec),'0'];
elseif length(num2str(timeDecSec)) == 4 %This if statments prevents errors in case time has 3 decimal places (Ex: 20.123us)
    fileName = [num2str(shot),'\',num2str(shot), ' ', locInput,' ',num2str(timeSec),'_',num2str(10*timeDecSec)];
else %Normal operation when time has 2 decimal places (Ex: 20.12us)
   fileName = [num2str(shot),'\',num2str(shot), ' ', locInput,' ',num2str(timeSec),'_',num2str(timeDecSec)];
end

% if ~exist(fileName,'var')
%     %throw(MException('MyComponent:noSuchVariable','Variable location not found :(');
%     return
% end

%File can include header or not, readtable() will handle non-numeric
%entries

% try
%     A = readtable(fileName);
% catch ME
%     disp('return')
%     visual = [0,0];
%     peakCenterMatrix=[0,0];
%     t = [0,0];
%     tplus1 = [0,0];
%     tminus1 = [0,0];
%     return
% end
% disp('did not return')

A = readtable(fileName);
A = A{:,:}; %This command converts to a matrix of doubles
[R,~] = size(A);
pixels = (1:R)';
%% 
%plot raw data from csv file
figure(1)
plot(pixels,A)

%Define window limits to exclude reflections from walls
windowLim1 = windowInput1;%86;<--limits for Aug 4 %189 <-- limits for Aug 11
windowLim2 = windowInput2; %224; %267
windowLim3 = windowInput3; %554; %495
windowLim4 = windowInput4; %710; %592

if shot == 210804034 || shot == 210804040 %this if statement may be unncessary since window inputs are now included in function arguments
    windowLim1 = 86;%<--limits for Aug 4 %189 <-- limits for Aug 11
    windowLim2 = 224; %267
    windowLim3 = 554; %495
    windowLim4 = 710;
end

%set all points under intensityCutoff to 0 to clean up plot
%intensityCutoff = input('Please enter intensity cutoff value: ');
%intensityCutoff = intensityInput;%880;

%%%%%%%%%Create code to replace user input for intensityCutoff, uncomment
%%%%%%%%%command above (input command) until code below works
conditionsMetCounter = 0;
intensityCutoff = 0;
for j = 950:-2:400
    %disp('New Loop')
    
    conditionsMetMatrix = zeros(1,1);
    for k = windowLim1:windowLim2
        %Condition: pixel must be below cutoff and next sequential pixel
        %must be above or pixel must be above and next one is below.
        %However, if both conditons are met (not just one), we will not
        %allow that because it means there is a single point as a peak.
                
        if (A(k)<=j && A(k+1) > j) || (A(k)>=j && A(k+1) < j) %&& ~((A(k)>=j && A(k+1)<j && A(k-1)<j) || (A(k)<=j && A(k+1)>j && A(k-1)>j)) %<-still working on this last part
            conditionsMetCounter = conditionsMetCounter + 1;
            conditionsMetMatrix(conditionsMetCounter) = k;
            %disp(k)
        end
    end
    
    for k = windowLim3:windowLim4
        %Condition: pixel must be below cutoff and next sequential pixel
        %must be above or pixel must be above and next one is below.
        %However, if both conditons are met (not just one), we will not
        %allow that because it means there is a single point as a peak.
                
        if (A(k)<=j && A(k+1) > j) || (A(k)>=j && A(k+1) < j) %&& ~((A(k)>=j && A(k+1)<j && A(k-1)<j) || (A(k)<=j && A(k+1)>j && A(k-1)>j)) %<-still working on this last part
            conditionsMetCounter = conditionsMetCounter + 1;
            conditionsMetMatrix(conditionsMetCounter) = k;
            %disp(k)
        end
    end
    
    
    if conditionsMetCounter == 4 && conditionsMetMatrix(1)<=windowLim2 && conditionsMetMatrix(2)<=windowLim2 && conditionsMetMatrix(3)>=windowLim3 && conditionsMetMatrix(4)>=windowLim3 %&& length(conditionsMetMatrix) == 4
        intensityCutoff = j;
        break
    else
        %disp(conditionsMetCounter);
        conditionsMetCounter = 0;
        clear conditionsMetMatrix
        %disp('New loop')
        %clc
    end
    
end
% for j = 950:-2:400
%     for k = 1:R-1
%         if (A(k) <= j && A(k+1) > j) || (A(k) >= j && A(k+1) < j)
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%These two for loops will 'clean up' the data only in the pixels that
%contain the viewing window
if intensityCutoff == 0 %This will handle any cases in which the above loops couldn't find a cutoff value
    intensityCutoff = input('Please enter intensity cutoff value: ');
end

intensityMatrix = zeros(R,1);
for k = windowLim1:windowLim2 %rows
    if A(k) >= intensityCutoff
        intensityMatrix(k) = A(k);
    end
end
    
for k = windowLim3:windowLim4 %rows
    if A(k) >= intensityCutoff
        intensityMatrix(k) = A(k);
    end
end

%% Temp section for 210804034 54.25us
%intensityMatrix(225:300) = 0;

%% Calculate pinch location/width
%plot cleaned up data
figure(2)
hold on
plot(pixels,intensityMatrix)
yline(intensityCutoff)
hold off

%determine peak width and location
counter = 0; %0 for left edge of peak, 1 for right edge of peak
peakCounter = 1; %will hold the value of how many peaks have been detected (should only be max of 2)
tempPeakValue = 0; %will hold the value for the pixel on the left edge of peak
%peakWidth = 0; %will hold the value of the width of the peak

    for k = 1:R %rows          
        if intensityMatrix(k) ~= 0
            
            if intensityMatrix(k-1)==0 && intensityMatrix(k+1)==0 %This should handle single pixel peaks
                peakCenterMatrix(peakCounter) = k;
                peakCounter = peakCounter + 1;
            end
            
            if (intensityMatrix(k-1)==0 || intensityMatrix(k+1)==0) && ~(intensityMatrix(k-1)==0 && intensityMatrix(k+1)==0) %This last part should disregard any single pixel peaks
                if counter == 0
                    tempPeakValue = k;
                    counter = 1;
                else
                    peakCenter = (tempPeakValue + k)/2;
                    %peakWidth = k-tempPeakValue;
                    counter = 0;
                    %peakWidthMatrix(peakCounter) = peakWidth;
                    peakCenterMatrix(peakCounter) = peakCenter;
                    peakCounter = peakCounter + 1;
                end
            end
            
        end               
    end
%     counter = 0;
%     peakCounter = 1;
%     tempPeakValue = 0;
%     peakWidth = 0;
%     peakCenter = 0;

%disp(peakWidthMatrix)
%disp('')
%disp(peakCenterMatrix)

%Calculate pinch location in chamber
upperWindowCenter = (windowLim2+windowLim1)/2;
upperWindowWidth = windowLim2-windowLim1;

lowerWindowCenter = (windowLim4+windowLim3)/2;
lowerWindowWidth = windowLim4-windowLim3;

%Plot cleaned up data with window limits shown
figure(3)
hold on
plot(pixels,intensityMatrix)
xline(upperWindowCenter,'LineWidth',1)
xline(windowLim1,'LineWidth',1)
xline(windowLim2,'LineWidth',1)

xline(lowerWindowCenter,'LineWidth',1)
xline(windowLim3,'LineWidth',1)
xline(windowLim4,'LineWidth',1)
hold off

windowWidth = 5; %centimeters
%positive X is left if looking downstream (towards hallway)
%a higher pixel indicates negative X direction for top view
%positive Y is up if looking downstream

if length(peakCenterMatrix) == 1 %Added this if else statement to handle times where I only need to examine the bottom window, this also assumes only the bottom window has a pinch instead of vice versa
    pinchXLocation = 3.0;
    pinchYLocation = (lowerWindowCenter - peakCenterMatrix(1))/lowerWindowWidth * windowWidth;
else
    pinchXLocation = (upperWindowCenter - peakCenterMatrix(1))/upperWindowWidth * windowWidth;
    pinchYLocation = (lowerWindowCenter - peakCenterMatrix(2))/lowerWindowWidth * windowWidth;
end

if abs(pinchXLocation) >= 2.5 || abs(pinchYLocation) >= 2.5
    disp(pinchXLocation)
    disp(pinchYLocation)
    disp(locInput)
    disp('\n')
end

%plot X and Y location
outerElectrodeRadius = 21.59/2; %cm %Double check these values
innerElectrodeRadius = 5.08; %cm

figure(4)
hold on

%Plot the pinch location along with the electrodes
th = 0:pi/50:2*pi;
xunit = outerElectrodeRadius * cos(th);
yunit = outerElectrodeRadius * sin(th);
xunit2 = innerElectrodeRadius * cos(th);
yunit2 = innerElectrodeRadius * sin(th);

plot(xunit, yunit);
plot(pinchXLocation,pinchYLocation,'x','MarkerSize',12,'LineWidth',2)

plot(xunit2, yunit2);
title(fileName)
grid on
hold off

%% Plot visual and jScope data together
%%%% Next step: correlating data w/ mode/probe location data
if compareToProbe == true
    %Connect to MDSplus and grab data
    mdsconnect('zappa.zap');
    mdsopen('zaphd',shot);
    MagTime = double(mdsvalue('dim_of(.5 * \R_OUTER * SIGMUL(\M_1_P15_NORM, COS(\PHI_M_1_P15 - $PI / 8)))'));
    MagTime = round(MagTime,10);
    
    %Probe values will be corrected for magnitude w/ the 100 in
    %front
    %XLoc_probes = -100*mdsvalue('.5 * \R_OUTER * SIGMUL(\M_1_P15_NORM, COS(\PHI_M_1_P15 - $PI / 8))');
    %YLoc_probes = -100*mdsvalue('.5 * \R_OUTER * SIGMUL(\M_1_P15_NORM, SIN(\PHI_M_1_P15 - $PI / 8))');
    XLoc_probes = 100*mdsvalue('\X_P15');
    YLoc_probes = 100*mdsvalue('\Y_P15');
    mdsclose;
    mdsdisconnect;


%jScope updates every 0.25us, so certain times from Kirana will not line up
%This if statement rounds down to the nearest 0.25us
if mod(time,0.25e-6) ~=0
    time = round(time - mod(time,0.25e-6),7);
end

correctedTime1 = round(time+1e-6,10);%time + 1us
correctedTime2 = round(time-1e-6,10);%time - 1us

%If using MDSplus, we only need to pull one time matrix
R1 = find(MagTime == time);
%R2 = find(yTime == time);

R3 = find(MagTime == correctedTime1);
%R4 = find(yTime == correctedTime1);

R5 = find(MagTime == correctedTime2);
%R6 = find(yTime == correctedTime2);

%Plot visual data, probe data @ t = -1,0,+1us, and electrodes
figure(5)
hold on
plot(xunit, yunit);
plot(xunit2, yunit2);

plot(pinchXLocation,pinchYLocation,'x','MarkerSize',12,'LineWidth',2)

plot(XLoc_probes(R1),YLoc_probes(R1),'x','MarkerSize',12,'LineWidth',2)%Plot location at t+0us
%plot(XLoc_probes(R3),YLoc_probes(R3),'x','MarkerSize',12,'LineWidth',2)%Plot location at t+1us
%plot(XLoc_probes(R5),YLoc_probes(R5),'x','MarkerSize',12,'LineWidth',2)%Plot location at t-1us

legend('Outer Electrode','Inner Electrode','Kirana Data','B Field Data');%,'time+1us','time-1us')
%title(strcat(num2str(shot),'---',num2str(timeSec),'.',num2str(timeDecSec),'us'))
title(strcat('Pulse',{' '},num2str(shot), ' at',{' '}, num2str(timeSec),'.',num2str(timeDecSec),'us (Upstream View)'))
xlabel('X Axis (cm)')
ylabel('Y Axis (cm)')
ax=gca;
ax.FontSize=16;
grid on
hold off

end
%% Determine phase and magnitude of locatons
%%%%Will need to adjust this code to remove R2,R4,R6,etc. if MDSplus works
% Phase
%0 degrees is at the top
%90 degrees is left (-X,0)
% phaseVisual = rad2deg(-atan2(pinchXLocation,pinchYLocation));
% phaseTime = rad2deg(atan2(100*XLoc_probes(R1,2),-100*YLoc_probes(R2,2)));
% phaseTimePlus1 = rad2deg(atan2(100*XLoc_probes(R3,2),-100*YLoc_probes(R4,2)));
% phaseTimeMinus1 = rad2deg(atan2(100*XLoc_probes(R5,2),-100*YLoc_probes(R6,2)));
% 
% magnitudeVisual = sqrt(pinchXLocation^2 + pinchYLocation^2);
% magnitudeTime = sqrt((100*XLoc_probes(R1,2))^2+(100*YLoc_probes(R2,2))^2);
% magnitudeTimePlus1 = sqrt((100*XLoc_probes(R3,2))^2+(100*YLoc_probes(R4,2))^2);
% magnitudeTimeMinus1 = sqrt((100*XLoc_probes(R5,2))^2+(100*YLoc_probes(R6,2))^2);
%% Determine Wavelength and plot locations against wavelength
%This task will be handled in wavelengthCalculator2()

% Take multiple y slices and curvefit them against a sine curve



% if wavelengthPos > 0
%     if wavelengthPos == 1
%         fileNameHoriz = [fileName,' Horiz Top'];
%         wavelength = wavelengthCalculator(fileNameHoriz,peakWidthMatrix(1),intensityInputHoriz); %1 for top window, 2 for bottom
%     else
%         fileNameHoriz = [fileName,' Horiz Bot'];
%         wavelength = wavelengthCalculator(fileNameHoriz,peakWidthMatrix(2),intensityInputHoriz); %1 for top window, 2 for bottom
%     end
% else
%     wavelength = 0;
% end
% %disp(wavelength);
% 
% % Convert from wavelength in pixels to cm
% if wavelengthPos == 1
%     wavelengthCM = wavelength/upperWindowWidth * 5;
% elseif wavelengthPos == 2
%     wavelengthCM = wavelength/lowerWindowWidth * 5;
% else
%     wavelengthCM = 0;
% end
% disp(wavelengthCM)
% 
%% Output variables
visual = [pinchXLocation, pinchYLocation]; %X & Y Location from visual data

if compareToProbe == false
    t = [0,0];
    tplus1 = [0,0];
    tminus1 = [0,0];
else
    t = [XLoc_probes(R1),YLoc_probes(R1)]; %X & Y Locaton from B probes @ same time as visual data
    tplus1 = [XLoc_probes(R3),YLoc_probes(R3)];% "   " @ t+1us from visual data
    tminus1 = [XLoc_probes(R5),YLoc_probes(R5)];% "  " @ t-1us from visual data 
end
end