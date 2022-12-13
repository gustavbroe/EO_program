%% Initializing MATLAB
clear; clc; close all;

% Adding function folder to path
addpath('functions', 'use')

% Loading materials information
MI = material_info(1);

% Varibale names for future tables
VarNames = {'f_ck','ApMoment','MomentCap','UtiFactor','UtiChange',...
    'GWP','PerChange'};

%----------------------- Span of strength classes -----------------------

% Range of strength classes and associated EPD
% (The letter indicates the cement type, R = RAPID)
f_ck = [16, 20, 25, 30, 35, 40]';
EPD_c = 'c' + string(f_ck) + 'R';

% Steel EPD
% (rp = Reinforcement products from Celsa Steel Service)
EPD_s = 'rp';


%------------------------------ Column SE01 ------------------------------

% Printing header
fprintf('COLUMN SE01:\n\n')

% Cross section dimensions
% (Information about the transverse reinforcement is found below)
h = 240;                % Height                        [mm]
w = 360;                % Width                         [mm]
                        % LONGITUDINAL
o_L = [16, 16, 16];     % Diameter                      [mm]
n_L = [5, 2, 5];        % Number of bars in each row    [-]
d_L = 145;              % Depth of middle of row        [mm]

% Allocation space in arrays
z = zeros(size(f_ck, 1), 1);
GWP = z;   ApMoment = z;   MomentCap = z;

for i = 1:size(f_ck, 1)

% Print current concrete class
fprintf('C%.f | ', f_ck(i));

% Preliminary calculations with inputs
M = def_c(f_ck(i), h, w, o_L, n_L, d_L);

% Using program functions                   % Saving results to table
M.e_0 = geo_imperfections(1, M);            
M = bending_moment(1, M);                   ApMoment(i) = M.M_0Ed;
M.varphi_ef = effective_creep(1, M);
check_soe(1, M);
M.M_Rd = moment_capacity(1, M);             MomentCap(i) = M.M_Rd;

% The calc_GWP does not account for two transverse hoops
[G, G_c, G_s] = GWP_SE01(M, MI, EPD_s, EPD_c(i));
                                            GWP(i) = G;

fprintf('\n')

end

% Utilization factor
UtiF = ApMoment./MomentCap;

% Saving results to table
SE01 = table(f_ck, ApMoment, MomentCap, UtiF, (UtiF-UtiF(end-1))/UtiF(end-1), ...
    GWP, (GWP-GWP(end-1))./GWP(end-1), 'VariableNames', VarNames);

% Transmission
fprintf('\n')


%------------------------------ Wall VES01 ------------------------------

% Printing header
fprintf('WALL VES01:\n\n')

% Cross section dimensions
% (Information about the transverse reinforcement is found below)
h = 220;                % Height                        [mm]
w = 3600;               % Width                         [mm]
                        % LONGITUDINAL
s_L = 150;              % Spacing                       [mm]
o_L = [6, 6];           % Diameter                      [mm]
n_L = [w/s_L, w/s_L];   % Number of bars in each row    [-]
d_L = [];               % Depth of middle of row        [mm]

% Allocation space in arrays
z = zeros(size(f_ck, 1), 1);
GWP = z;   ApMoment = z;   MomentCap = z;

for i = 1:size(f_ck, 1)

% Print current concrete class
fprintf('C%.f | ', f_ck(i));

% Preliminary calculations with inputs
M = def_w(f_ck(i), h, w, o_L, n_L, d_L);

% Using program functions                   % Saving results to table
M.e_0 = geo_imperfections(1, M);            
M = bending_moment(1, M);                   ApMoment(i) = M.M_0Ed;
M.varphi_ef = effective_creep(1, M);
check_soe(1, M);
M.M_Rd = moment_capacity(1, M);             MomentCap(i) = M.M_Rd;
G = calc_GWP(2, [EPD_c(i), EPD_s], M);      GWP(i) = G;

fprintf('\n')

end

% Utilization factor
UtiF = ApMoment./MomentCap;

% Saving results to table
VES01 = table(f_ck, ApMoment, MomentCap, UtiF, (UtiF-UtiF(end-1))/UtiF(end-1), ...
    GWP, (GWP-GWP(end-1))./GWP(end-1), 'VariableNames', VarNames);

% Transmission
fprintf('\n')



% ----------------------- Printing/Saving results -----------------------

fprintf('Results for column SE01: \n');   disp(SE01)
fprintf('Results for wall VES01: \n');   disp(VES01)

% Saving to excel
writetable(SE01,'CementAnalysis.xlsx','Sheet',1);
writetable(VES01,'CementAnalysis.xlsx','Sheet',1,'Range','A10');

% -------------------------- GWP-GHG for column --------------------------

function [total, concrete, steel] = GWP_SE01(M, MI, EPD_s, EPD_c)

% Transverse reinforcement dimensions
o_T = 6;                % Diamter                       [mm]
w_T = 200;              % Width of middle hoop          [mm]
s_T = 150;              % Spacing between group of hoops[mm]

% Longitudinal
    % Volumen
    V_L = M.A_s * M.L;                                  %[m^3]

% Transverse
% (Assuming that the hoops are square)
    % Area of one transverse bars
    A_1T = pi * (o_T*1E-3/2)^2;                         %[m^2]
    % Length of hoops per set
    l_T = 4*(M.h-2*M.c) + 2*(M.w-2*M.c) + 2*w_T*1E-3;   %[m]
    % Number of hoops groups
    n_T = M.L/(s_T*1E-3);                               %[-]
    % Volume
    V_T = A_1T * l_T * n_T;                             %[m^3]

% GWP-GHG of steel
steel = (V_L + V_T)*MI.(EPD_s).GWP;                     %[kg CO2-eq]

% Concrete
    % Volumen
    V_c = M.w*M.h*M.L - (V_T + V_L);                    %[m^3]
    % GWP-GHG
    concrete = V_c*MI.(EPD_c).GWP;                      %[kg CO2-eq]

% GWP-GHG of SE01
total = concrete + steel;                               %[kg CO2-eq]

% Printing result
fprintf('GWP_(tot, s, c) = (%.2f, %.2f, %.2f) [kg CO2-eq.] | ', total, steel, concrete);

end