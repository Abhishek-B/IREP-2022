%{ 
IEEE 13 node impedence matrices
January 2022 - Abhishek B.
%}

%% imperial to metric unit conversion

% miles to KM
mi2km = 1.609344;

%% Config 601

R_601 = [0.3465, 0.1560, 0.1580;
     0.1560, 0.3375, 0.1535;
     0.1580, 0.1535, 0.3413];
 
X_601 = [1.0179, 0.5017, 0.4236;
     0.5017, 1.0478, 0.3849;
     0.4236, 0.3849, 1.0348];

Z_601 = (R_601+1i*X_601)*mi2km;

%% Config 602

R_602 = [0.7526, 0.1580, 0.1560;
         0.1580, 0.7475, 0.1535;
         0.1560, 0.1535, 0.7436];

X_602 = [1.1814, 0.4236, 0.5017;
         0.4236, 1.1983, 0.3849;
         0.5017, 0.3849, 1.2112];

Z_602 = (R_602+1i*X_602)*mi2km;

%% Config 603

R_603 = [0, 0, 0;
         0, 1.3294, 0.2066;
         0, 0.2066, 1.3238];

X_603 = [0, 0, 0;
         0, 1.3471, 0.4591;
         0, 0.4591, 1.3569];

Z_603 = (R_603+1i*X_603)*mi2km;

%% Config 604

R_604 = [1.3238, 0, 0.2066;
         0, 0, 0;
         0.2066, 0, 1.3294];

X_604 = [1.3569, 0, 0.4591;
         0, 0, 0;
         0.4591, 0, 1.3471];
     
Z_604 = (R_604+1i*X_604)*mi2km;

%% Config 605

R_605 = [0, 0, 0;
         0, 0, 0;
         0, 0, 1.3292];

X_605 = [0, 0, 0;
         0, 0, 0;
         0, 0, 1.3475];
     
Z_605 = (R_605+1i*X_605)*mi2km;

%% Config 606

R_606 = [0.7982, 0.3192, 0.2849;
         0.3192, 0.7891, 0.3192;
         0.2849, 0.3192, 0.7982];
     
X_606 = [0.4463, 0.0328, -0.0143;
         0.0328, 0.4041, 0.0328;
         -0.0143, 0.0328, 0.4463];
     
Z_606 = (R_606+1i*X_606)*mi2km;

%% Config 607

R_607 = [1.3425, 0, 0;
         0, 0, 0;
         0, 0, 0];
     
X_607 = [0.5124, 0, 0;
         0, 0, 0;
         0, 0, 0];
     
Z_607 = (R_607+1i*X_607)*mi2km;

%% Matrix R/X for example
% These are 39x39 matrices...but automating the task might take longer than
% just doing it manually since it is only for one example



