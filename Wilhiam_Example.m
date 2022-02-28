%% Example worked with Wilhiam

v0 = 2.4e3; % Voltage in Volts

t = 2; % time steps

p1 = [3e3;       % IN WATTS
     2e3;
     1e3];
 
q1 = [0;0;0];   % In VARS

p2 = [2e3;
      1e3;
      -6e3];
  
q2 = 0;

%%

% miles to feet
mi2ft = 1/5280;

% Config 601

R_601 = 2000*mi2ft*[0.3465, 0.1560, 0.1580;
     0.1560, 0.3375, 0.1535;
     0.1580, 0.1535, 0.3413];

%% Baseline votages

U_b1 = (v0^2)*ones(3,1) - R_601*p1;
U_b2 = (v0^2)*ones(3,1) - R_601*p2;

%%

figure()
plot(sqrt(U_b1)/v0, 'r')
hold on
grid on
plot(sqrt(U_b2)/v0, 'b')
hold on 
grid on
plot(ones(3,1), 'k')
ylim([0.9999,1.00006])
