clear all
close all
clc


shotInput = 211209004;
timeInput = 73.5E-6;
locInput = 'P15';
locInputPixel = 331;
windowInput1 = 219;
windowInput2 = 331;
windowInput3 = 562;
windowInput4 = 727;
compareToProbe = true;
%Analyze visual data
[visual, t, tplus1, tminus1,peakCenterMatrix] = Pinch_Location_Calculator(shotInput,timeInput,locInput,windowInput1,windowInput2,windowInput3,windowInput4,compareToProbe);
%Going to run Pinch_Location_Calculator again for 211118018 for the other
%time step, comment the next 2 lines out for all other pulses
% timeInput2 = 59.75E-6;
% [visual2, t, tplus1, tminus1,peakCenterMatrix] = Pinch_Location_Calculator(shotInput,timeInput2,locInput,windowInput1,windowInput2,windowInput3,windowInput4,compareToProbe);

%Pull B Field data
mdsconnect('zappa.zap');
mdsopen('zaphd',shotInput);

BField = double(100*mdsvalue('\Y_P15'));
MagTime = double(1E6*mdsvalue('dim_of(.5 * \R_OUTER * SIGMUL(\M_1_P15_NORM, COS(\PHI_M_1_P15 - $PI / 8)))'));
MagTime = round(MagTime,10);
mdsclose;
mdsdisconnect;

%Analyze interferometry data
close all
[YLocation,TimeArray,LocArray] = Interferometry_Pinch_Location2(shotInput,timeInput);

disp(['B Field Data: ',num2str(t(2))]) %jScope data
disp(['Visual Data: ',num2str(visual(2))]) %Visual data
disp(['IF Data: ',num2str(YLocation)]) %IF data

close all
figure(1)
hold on
plot(TimeArray,LocArray,'LineWidth',2)
%xline(timeInput*1e6,'r','LineWidth',2.5)
scatter(timeInput*1E6,visual(2),50,'filled','r')
%Comment out next line for any pulses except 211118018
% scatter(timeInput2*1E6,visual2(2),50,'filled','r')

legend('IF Centroid Location','Kirana Centroid Location')
xlabel('Time (us)')
ylabel('Centroid location (cm)')
title(['Pulse ',num2str(shotInput),' centroid location from IF and Kirana'])
ax=gca;
ax.FontSize=16;
hold off

figure(2)
hold on
plot(TimeArray,LocArray,'k','LineWidth',2)
scatter(timeInput*1E6,visual(2),50,'filled','r')
plot(MagTime,BField,'b','LineWidth',2)
plot(MagTime,BField+1,'c')
plot(MagTime,BField-1,'c')
%Comment out next line for any pulses except 211118018
% scatter(timeInput2*1E6,visual2(2),50,'filled','r')

legend('IF Centroid Location','Kirana Centroid Location','B Field Centroid Location','1 cm confidence margin for B Field')
xlabel('Time (us)')
ylabel('Centroid location (cm)')
title(['Pulse ',num2str(shotInput),' centroid location from IF, B Field Probes, and Kirana'])
ax=gca;
ax.FontSize=16;
axis([0,100,-4,4])
hold off
% figure(4)
% t = linspace(-1.5,1.5,100);
% x = timeInput*1e6;
% y = t;
% z = 1e21;
% line(x,y,z)
%9 Dec Kirana Window Settings
%P15 location: X = 331
%Top Window: Not Necessary
%Bottom Window: 562-727