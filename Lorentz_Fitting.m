% Lorentz fitting fcn
function [gamma,x0, fudgeFactor,RMSval]=Lorentz_Fitting(shot,timeInput,CH1,CH2,CH3,t)
%% Pull MDSplus data

%shot = 220204018;
% mdsconnect('zappa.zap');
% mdsopen('zaphd',shot);
% 
% CH1 = abs(mdsvalue('\NE_1'))'; %Top laser
% CH2 = abs(mdsvalue('\NE_2'))'; %Middle laser
% CH3 = abs(mdsvalue('\NE_3'))'./cos(deg2rad(5)); %Bottom laser 
% t = round(double(1e6*mdsvalue('dim_of(\NE_1)')'),3); %Convert from usec to sec make it easier
% 
% mdsclose;
% mdsdisconnect;

%% Find data at time timeInput
timeInput = timeInput*1e6;
if mod(timeInput,.04) ~= 0
    timeInput = timeInput - mod(timeInput,.04);
end
returnIndex = find(t==round(timeInput,2));
%% Fit Lorentz curve with fminsearch

%This is my Lorentz distribution function, I am trying to constrain X(2)
%such that -1.5<=X(2)=1.5, or abs(X(2))-1.5<=0.
y = @(X,loc)X(3)*1/pi*X(1)/((loc-X(2))^2+X(1)^2);

%These are single values, typically on the order of 1E20 to 1E23
Ch1 = double(CH1(returnIndex));
Ch2 = double(CH2(returnIndex));
Ch3 = double(CH3(returnIndex));
Ch = [Ch1,Ch2,Ch3];

%This is my objective function, the root mean squared equation.
%Essentially, I am trying to minimize this by adjusting the values of X in
%our Lorentz distribution.
RMS = @(X)sqrt(((y(X,1.5)-Ch1)^2+(y(X,0)-Ch2)^2 + (y(X,-1.5)-Ch3)^2)/3)+1E30*heaviside(abs(X(2))-1.5);

%This is my initial guess for the minimizing function
initial = [1,0,10^(floor(log10(max(Ch))))];%[gamma,x0,fudgeFactor]

%Trying to set different options to help get a better result
v=[.1,.01,1E20];
%options = optimset('PlotFcns',@optimplotfval,'MaxIter',600,'MaxFunEvals',1000,'FinDiffRelStep',v);%This shows the minimization over each iteration,not necessary for final code
options = optimset('MaxIter',600,'MaxFunEvals',1800,'Display','none');

[Output,RMSval] = fminsearch(RMS,initial,options);%Minimizing root mean squared eqn for Cauchy-Lorentz distribution fcn

if RMSval>1E17 %Rerunning if it didn't optimize completely on the first time
    initial=Output;
    [Output,RMSval] = fminsearch(RMS,initial,options);
end

gamma = Output(1);
x0 = Output(2);
fudgeFactor = Output(3);
%% Fit Lorentz curve with fmincon
for d=1
% y = @(X,loc)X(3)*1/pi*X(1)/((loc-X(2))^2+X(1)^2);
% Ch1 = double(CH1(returnIndex));
% Ch2 = double(CH2(returnIndex));
% Ch3 = double(CH3(returnIndex));
% Ch = [Ch1,Ch2,Ch3];
% 
% RMS = @(X)sqrt(((y(X,1.5)-Ch1)^2+(y(X,0)-Ch2)^2 + (y(X,-1.5)-Ch3)^2)/3);
% initial = [1,0,10^(floor(log10(max(Ch))))];%[gamma,x0,fudgeFactor]
% 
% function [c,ceq] = conditions(X)
% c = abs(X(2))-1.5; %|x0|<=1.5
% ceq=[];
% end
% 
% nonlcon = @conditions;
% v=[.1,.01,1E20];
% options = optimoptions(@fminimax,'PlotFcns',@optimplotfval,'MaxIter',600,'MaxFunEvals',1800,'FiniteDifferenceStepSize',v,'FunctionTolerance',1E10);
% %optimset('PlotFcns',@optimplotfval,'MaxIter',600,'MaxFunEvals',1800,'FinDiffRelStep',v);
% [Output] = fminimax(RMS,initial,[],[],[],[],[],[],nonlcon,options);
% gamma = Output(1);
% x0 = Output(2);
% fudgeFactor = Output(3);
end
%% Plot
% if abs(timeInput-50)<1E-5
%     
% close all
% locSpace = linspace(-1.5,1.5,100);
% hold on
% plot(locSpace,fudgeFactor*1/pi*gamma./((locSpace-x0).^2+gamma^2))
% scatter([1.5,0,-1.5],[Ch1,Ch2,Ch3])
% xline(x0)
% title(['Pulse ',num2str(shot), ' at ', num2str(timeInput),'us with Lorentz Dsitribution Curve Fit'])
% xlabel('Location (cm)')
% ylabel('Electron Number Density (m^{-3})')
% legend('Lorentz Distribution','Density Data Points')
% hold off
% ax=gca;
% ax.FontSize=16;
% end

end