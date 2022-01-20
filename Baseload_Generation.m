%{
Generation of Base Load. Files from Liz. Will need to clean up the process
eventually.
-January 2022, Abhishek B.
%}

clear all
close all

load('AusgridDataset_3yearv4.mat')

%% Edits Abby
% trying this now with aggregate customers and day 217 which should be the
% hottest day of that year

% Edit - Abby: trying this for day 217, which should be the peak day

gtot = zeros(1,48);
ltot = zeros(1,48);

p = zeros(70,48);
q = zeros(70,48);

theta = [];
theta(1) = acos(160/110);



for i=1:70 % using 70 customers
    cust = totalCustData(i); % customer i data
    custiday217 = cust.DayData(217);
    
    gtot = gtot + custiday217.GG*2;
    ltot = ltot + custiday217.GC*2;
    
    p217 = custiday217.GC*2 - custiday217.GG*2;
    
    q217 = (1/cos(theta))*p217;
    
    p(i,:) = p217;
    q(i,:) = q217;
end

ptot = zeros(1,48);
qtot = zeros(1,48);

for i=1:70
    ptot = ptot + p(i,:);
    qtot = qtot + q(i,:);
end

figure()
plot(ptot)
hold on
plot(qtot)
grid on
legend('p','q')

%% Creating voltage profile from p and q.

load('RX.mat')
