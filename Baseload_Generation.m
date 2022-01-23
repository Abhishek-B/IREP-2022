%{
Generation of Base Load. Files from Liz. Will need to clean up the process
eventually.
-January 2022, Abhishek B.
%}

clear all
close all

load('AusgridDataset_3yearv4.mat')

%%

day = 217;

%% Nodes 1 & 2
% No spot load data is given, so I assume the base load here to be zero

p1a = zeros(1,48);
q1a = zeros(1,48);
p1b = zeros(1,48);
q1b = zeros(1,48);
p1c = zeros(1,48);
q1c = zeros(1,48);
p2a = zeros(1,48);
q2a = zeros(1,48);
p2b = zeros(1,48);
q2b = zeros(1,48);
p2c = zeros(1,48);
q2c = zeros(1,48);

%% Node 3a
p3a = zeros(1,48);
theta3a = acos(160/110);

for i=1:70 % using 70 customers
    cust = totalCustData(i); % customer i data
    custiday217 = cust.DayData(day);
    p3a = p3a + (custiday217.GC*2 - custiday217.GG*2);
end
q3a = (1/cos(theta3a))*p3a;

%% Node 3b
p3b = zeros(1,48);
theta3b = acos(120/90);

for i=1:52
    cust = totalCustData(i); % customer i data
    custiday217 = cust.DayData(day);
    p3b = p3b + (custiday217.GC*2 - custiday217.GG*2);
end
q3b = (1/cos(theta3b))*p3b;

%% Node 3c
p3c = zeros(1,48);
theta3c = acos(120/90);

for i=1:52
    cust = totalCustData(i); % customer i data
    custiday217 = cust.DayData(day);
    p3c = p3c + (custiday217.GC*2 - custiday217.GG*2);
end
q3c = (1/cos(theta3c))*p3c;

%% Node 4b
% Node 4 phase A, the IEEE stuff has spot load as 0...not sure whats 
% happening there

p4b = zeros(1,48);
theta4b = acos(170/125);

for i=1:73
    cust = totalCustData(i); % customer i data
    custiday217 = cust.DayData(day);
    p4b = p4b + (custiday217.GC*2 - custiday217.GG*2);
end
q4b = (1/cos(theta4b))*p4b;

%% Node 4c
p4c = zeros(1,48);
q4c = zeros(1,48);

%% Node 5b
% Node 5 phase A and C, the IEEE stuff has spot load as 0...not sure whats 
% happening there

p5b = zeros(1,48);
theta5b = acos(230/132);

for i=1:100
    cust = totalCustData(i); % customer i data
    custiday217 = cust.DayData(day);
    p5b = p5b + (custiday217.GC*2 - custiday217.GG*2);
end
q5b = (1/cos(theta5b))*p5b;

%% Node 5c
p5c = zeros(1,48);
q5c = zeros(1,48);

%% Node 6a

p6a = zeros(1,48);
theta6a = acos(385/220);

for i=1:175
    cust = totalCustData(i); % customer i data
    custiday217 = cust.DayData(day);
    p6a = p6a + (custiday217.GC*2 - custiday217.GG*2);
end
q6a = (1/cos(theta6a))*p6a;

%% Node 6b and 6c
p6b = p6a;
q6b = q6a;
p6c = p6a;
q6c = q6a;

%% Node 7a & 7b & 7c

p7a = zeros(1,48);
q7a = zeros(1,48);
p7b = zeros(1,48);
q7b = zeros(1,48);

p7c = zeros(1,48);
theta7c = acos(170/151);

for i=1:73
    cust = totalCustData(i); % customer i data
    custiday217 = cust.DayData(day);
    p7c = p7c + (custiday217.GC*2 - custiday217.GG*2);
end
q7c = (1/cos(theta7c))*p7c;

%% Node 8a

p8a = zeros(1,48);
theta8a = acos(485/190);

for i=1:227
    cust = totalCustData(i); % customer i data
    custiday217 = cust.DayData(day);
    p8a = p8a + (custiday217.GC*2 - custiday217.GG*2);
end
q8a = (1/cos(theta8a))*p8a;

%% Node 8b

p8b = zeros(1,48);
theta8b = acos(68/60);

for i=1:32
    cust = totalCustData(i); % customer i data
    custiday217 = cust.DayData(day);
    p8b = p8b + (custiday217.GC*2 - custiday217.GG*2);
end
q8b = (1/cos(theta8b))*p8b;

%% Node 8c

p8c = zeros(1,48);
theta8c = acos(290/212);

for i=1:133
    cust = totalCustData(i); % customer i data
    custiday217 = cust.DayData(day);
    p8c = p8c + (custiday217.GC*2 - custiday217.GG*2);
end
q8c = (1/cos(theta8c))*p8c;

%% Node 9a & 9c

p9a = zeros(1,48);
q9a = zeros(1,48);
p9c = zeros(1,48);
q9c = zeros(1,48);

%% Node 10c

p10c = zeros(1,48);
theta10c = acos(170/80);

for i=1:73
    cust = totalCustData(i); % customer i data
    custiday217 = cust.DayData(day);
    p10c = p10c + (custiday217.GC*2 - custiday217.GG*2);
end
q10c = (1/cos(theta10c))*p10c;

% figure()
% plot(p10c)
% hold on
% plot(q10c)
% grid on
% legend('p','q')

%% Node 11a

p11a = zeros(1,48);
theta11a = acos(128/86);

for i=1:56
    cust = totalCustData(i); % customer i data
    custiday217 = cust.DayData(day);
    p11a = p11a + (custiday217.GC*2 - custiday217.GG*2);
end
q11a = (1/cos(theta11a))*p11a;

%% Node 12a,b,c
p12a = zeros(1,48);
q12a = zeros(1,48);
p12b = zeros(1,48);
q12b = zeros(1,48);
p12c = zeros(1,48);
q12c = zeros(1,48);


%% Final p and q vectors (matrices really) over the time period

P = [ p1a ; p1b ; p1c ; p2a ; p2b ; p2c ; p3a ; p3b ; p3c ; p4b ; p4c ; ...
     p5b ; p5c ; p6a ; p6b ; p6c ; p7a ; p7b ; p7c ; p8a ; p8b ; p8c ; ...
     p9a ; p9c ; p10c ; p11a ; p12a ; p12b ; p12c];

Q = [ q1a ; q1b ; q1c ; q2a ; q2b ; q2c ; q3a ; q3b ; q3c ; q4b ; q4c ; ...
     q5b ; q5c ; q6a ; q6b ; q6c ; q7a ; q7b ; q7c ; q8a ; q8b ; q8c ; ...
     q9a ; q9c ; q10c ; q11a ; q12a ; q12b ; q12c];


%% Creating voltage profile from p and q.

v0 = 4.16*ones(29,48);

load('RX.mat')

V = v0 - R*P - X*Q;

V_base = (1/4.16)*V;

% Saving the base load
save('baseload.mat', 'V_base');

