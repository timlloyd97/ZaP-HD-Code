function [YLocation,TimeArray,LocArray] = Interferometry_Pinch_Location2(shot,timeInput)
%% Pull data from MDSplus and tailor it from 1 to 100us
mdsconnect('zappa.zap');
mdsopen('zaphd',shot);

CH1 = abs(mdsvalue('\NE_1'))'; %Top laser
CH2 = abs(mdsvalue('\NE_2'))'; %Middle laser
CH3 = abs(mdsvalue('\NE_3'))'./cos(deg2rad(5)); %Bottom laser 

% CH1 = (mdsvalue('\NE_1'))'; %Top laser
% CH2 = (mdsvalue('\NE_2'))'; %Middle laser
% CH3 = (mdsvalue('\NE_3'))'./cos(deg2rad(5)); %Bottom laser 

%Going to just try adjusting everything vertically so they all start at 0
% CH1 = abs(CH1-CH1(1));
% CH2 = abs(CH2-CH2(1));
% CH3 = abs(CH3-CH3(1));

BProbes = 100*mdsvalue('\Y_P15');
%BProbes = 100*mdsvalue('.5 * \R_OUTER * SIGMUL(\M_1_P15_NORM, SIN(\PHI_M_1_P15 - $PI / 8))');
t = round(double(1e6*mdsvalue('dim_of(\NE_1)')'),2); %Convert from us to make it easier
t_BProbes = round(1e6*mdsvalue('dim_of(\Y_P15)'),5);

mdsclose;
mdsdisconnect;

t0 = find(t==0);
t100 = find(t==100);

t0_B = find(t_BProbes == 0);
t100_B = find(t_BProbes == 100);

%Reassigning data to only show between 0 and 100us
CH1 = CH1(t0:t100);
CH2 = CH2(t0:t100);
CH3 = CH3(t0:t100);
%BProbes = BProbes(t0:t100);
t = t(t0:t100);

BProbes = BProbes(t0_B:t100_B);
t_BProbes = t_BProbes(t0_B:t100_B);

CH1loc = 1.5;
CH2loc = 0;
CH3loc = -1.5;

t_test = double([t,t,t]);
Ne = double([CH1,CH2,CH3]);
loc = [CH1loc*ones(1,length(t)),CH2loc*ones(1,length(t)),CH3loc*ones(1,length(t))];

%[xq,yq] = meshgrid(double(t(t0)):.04:double(t(t100)),-1.5:.5:1.5);
%^This line works, I'm just experimenting with how many values I can plot
[xq,yq] = meshgrid(double(t(1)):.04:double(t(end)),-1.5:.5:1.5);
vq = griddata(t_test,loc,Ne,xq,yq);
figure(1)
mesh(xq,yq,vq)
xlabel('Time (us)')
ylabel('Location (cm)')
zlabel('Density (m^{-3})')
title(['Pulse ',num2str(shot),' density surface plot'])
ax=gca;
ax.FontSize=16;
%% Begin fitting the data to Lorentz distribution
x0 = zeros(2501,1);
gamma = zeros(2501,1);
fudgeFactor = zeros(2501,1);

count=1;
for j = .0:.04:100
    [gamma(count),x0(count), fudgeFactor(count),~]=Lorentz_Fitting(shot,j*1E-6,CH1,CH2,CH3,t);
    count=count+1;
end

if mod(timeInput*1E6,0.04) ~=0
    timeInput = round(timeInput*1E6 - mod(timeInput*1E6,0.04),7);
else
    timeInput=round(timeInput*1E6,7);
end
returnIndex = find(t==timeInput);
TimeArray = t';
YLocation = x0(returnIndex);%x0(returnIndex); %Need to fix this to handle decimal values


%Will adjust LocArray vertically based on average of first 10 values
offset = mean(x0(1:10));
if abs(offset) <=1.5
    LocArrayOffset = x0-offset;
else
    LocArrayOffset = x0;
end
LocArray = x0;

figure(2)
plot(0:.04:100,LocArray,t_BProbes,BProbes,'LineWidth',2)%0:.04:100,LocArrayOffset,
legend('IF Centroid Location','B Field Centroid Location')%'Adjusted IF Centroid Location',
xlabel('Time (us)')
ylabel('Centroid location (cm)')
title(['Pulse ', num2str(shot), ' centroid location from IF and B Field Data'])
ax2=gca;
ax2.FontSize=16;

%% Plot centroid data on surface plot
Ne_fixed = [CH1;CH2;CH3]; %Same Ne matrix, but CH1-3 are vertically stacked as opposed to horizontally concatenated
[Ne_max] = max(abs(Ne_fixed));%These matrices hold the max value at every time step as well as which channel they're on
figure(3)
hold on
%plot3(t,LocArrayOffset,Ne_max,'LineWidth',2);
plot3(t,LocArray,Ne_max,'LineWidth',2,'Color','r');
mesh(xq,yq,vq)
xlabel('Time (us)')
ylabel('Location (cm)')
zlabel('Density (m^{-3})')
legend('IF Centroid Location')%'Adjusted IF Centroid Location',
title(['Pulse ',num2str(shot),' density surface plot with IF Centroid Data'])
ax3=gca;
ax3.FontSize=16;
hold off