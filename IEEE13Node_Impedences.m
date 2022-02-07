%{ 
IEEE 13 node impedence matrices
January 2022 - Abhishek B.
%}

%% imperial to metric unit conversion

% miles to feet
mi2ft = 1/5280;

%% Config 601

R_601 = [0.3465, 0.1560, 0.1580;
     0.1560, 0.3375, 0.1535;
     0.1580, 0.1535, 0.3413];
 
X_601 = [1.0179, 0.5017, 0.4236;
     0.5017, 1.0478, 0.3849;
     0.4236, 0.3849, 1.0348];

Z_601 = (R_601+1i*X_601)*mi2ft;

%% Config 602

R_602 = [0.7526, 0.1580, 0.1560;
         0.1580, 0.7475, 0.1535;
         0.1560, 0.1535, 0.7436];

X_602 = [1.1814, 0.4236, 0.5017;
         0.4236, 1.1983, 0.3849;
         0.5017, 0.3849, 1.2112];

Z_602 = (R_602+1i*X_602)*mi2ft;

%% Config 603

R_603 = [0, 0, 0;
         0, 1.3294, 0.2066;
         0, 0.2066, 1.3238];

X_603 = [0, 0, 0;
         0, 1.3471, 0.4591;
         0, 0.4591, 1.3569];

Z_603 = (R_603+1i*X_603)*mi2ft;

%% Config 604

R_604 = [1.3238, 0, 0.2066;
         0, 0, 0;
         0.2066, 0, 1.3294];

X_604 = [1.3569, 0, 0.4591;
         0, 0, 0;
         0.4591, 0, 1.3471];
     
Z_604 = (R_604+1i*X_604)*mi2ft;

%% Config 605

R_605 = [0, 0, 0;
         0, 0, 0;
         0, 0, 1.3292];

X_605 = [0, 0, 0;
         0, 0, 0;
         0, 0, 1.3475];
     
Z_605 = (R_605+1i*X_605)*mi2ft;

%% Config 606

R_606 = [0.7982, 0.3192, 0.2849;
         0.3192, 0.7891, 0.3192;
         0.2849, 0.3192, 0.7982];
     
X_606 = [0.4463, 0.0328, -0.0143;
         0.0328, 0.4041, 0.0328;
         -0.0143, 0.0328, 0.4463];
     
Z_606 = (R_606+1i*X_606)*mi2ft;

%% Config 607

R_607 = [1.3425, 0, 0;
         0, 0, 0;
         0, 0, 0];
     
X_607 = [0.5124, 0, 0;
         0, 0, 0;
         0, 0, 0];
     
Z_607 = (R_607+1i*X_607)*mi2ft;

%% Matrix R/X for example
% These are 29x29 matrices...but automating the task might take longer than
% just doing it manually since it is only for one example

% forget what I said before, I must've been high...841 manual computations
% are ridiculous 

% Phases
a=1;
b=2;
c=3;


% Indices of the matrices
index = {[1,a], [1,b], [1,c], [2,a], [2,b], [2,c], [3,a], [3,b], [3,c],...
         [4,b], [4,c], [5,b], [5,c], [6,a], [6,b], [6,c], [7,a], [7,b],...
         [7,c], [8,a], [8,b], [8,c], [9,a], [9,c], [10,c], [11,a],...
         [12,a], [12,b], [12,c]};

% Edges along nodes

E={};
E{1} = ["01"];
E{2} = ["01", "12"];
E{3} = ["01", "12", "23"];
E{4} = ["01", "14"];
E{5} = ["01", "14", "45"];
E{6} = ["01", "16"];
E{7} = ["01", "16", "67"];
E{8} = ["01", "16", "67", "78"];
E{9} = ["01", "16", "69"];
E{10} = ["01", "16", "69", "910"];
E{11} = ["01", "16", "69", "911"];
E{12} = ["01", "16", "612"];

% Dictionary to translate edge to IEEE configs

% These are ignored in the EV models, so I will assign them a different 60X
% config from IEEE 13 node
transformer=Z_602; % same as edge 1-2
Swtch=Z_606; % same as edge 7-8
% decided to replace them with the impedences of the adjacent lines

dict = containers.Map();

% Here I multiply the impedence given with the length of the particular
% line also give. This seems important since some lines have different
% length, but same config for impedence.

dict("01") = Z_601*2000;
dict("12") = Z_602*500;
dict("23") = transformer*500; % note IEEE has this lenght as 0, but they
                              % also have it as a transformer
dict("14") = Z_603*500;
dict("45") = Z_603*300;
dict("16") = Z_601*2000;
dict("67") = Swtch*500; % note IEEE has this lenght as 0, but they
                        % also have it as a switch
dict("78") = Z_606*500;
dict("69") = Z_604*300;
dict("910") = Z_605*300;
dict("911") = Z_607*800;
dict("612") = Z_601*1000;


% building R/X

R = zeros(29,29);
X = zeros(29,29);
omega = exp(-2*pi*1i/3);

for i=1:length(index) %row entry
    for j=1:length(index) %col entry
        
        row = index{i};
        col = index{j};
        
        edgei = E{row(1)};
        edgej = E{col(1)};
        
        common_edges = intersect(edgei, edgej);
        
        Impedence_sum = zeros(3,3);
        for k=1:length(common_edges)
            Impedence_sum = Impedence_sum + dict(common_edges(k));
        end
        
        entry = Impedence_sum(row(2),col(2));
        coeff = omega^(row(2)-col(2));
        
        R(i,j) = 2*real(conj(entry)*coeff);
        X(i,j) = -2*imag(conj(entry)*coeff);
    end
end


% Saving the matrices
save('RX.mat', 'R', 'X');