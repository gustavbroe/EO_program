%% Initializing Matlab
clear, clc, close all hidden;

% Adding function folder to path
addpath('functions')

% Width of wall
w = 3600;           %[mm]

% Reinforcement steel EPD
EPD_s = 'rp';

% Starting counter
count = 1;

%% Control scenario
fprintf('\n---------- Control ---------- \n');

% Characteristic compressive strength of concrete
% (16, 20, 25, 30 [F/R]  |  35, 40 [L/R])
f_ck = 35;   EPD_c = 'c35R';

% Wall thickness
h = 220;

% Diameter of longitunal reinforcement
o_L = 6;

% Spacing between longitunal reinforcement
s_L = 150;

% Printing header
fprintf('\n(C%.f/h%.f/ø%.f/s%.f):   ', f_ck, h, o_L, s_L)

% Calculations
[M_Ed, M_Rd, GWP] = mom_dem(f_ck, h, w, [o_L, o_L], [w/s_L, w/s_L], [], EPD_c, EPD_s);

% Saving results
res(count,:) = [f_ck, h, o_L, s_L, GWP, M_Rd, M_Ed.sm, M_Ed.nc, M_Ed.ns];
count = count + 1;



%% Scenario A
fprintf('\n\n\n---------- Scenario A ---------- \n');

% Printing header
fprintf('\nChange of concrete strength class: \n');

% Characteristic compressive strength of concrete
% (16, 20, and 35 is currently available)
f_ck = [35 20 16];   EPD_c = {'c35R', 'c20F', 'c16F'};

% Wall thickness
h = 220;

% Diameter of longitunal reinforcement
o_L = 6;

% Spacing between longitunal reinforcement
s_L = 150;

% Looping through each case
for i = 1:size(f_ck, 2)

    % Printing header
    fprintf('\n(C%.f/h%.f/ø%.f/s%.f):   ', f_ck(i), h, o_L, s_L)
    
    % Calculations
    [M_Ed, M_Rd, GWP] = mom_dem(f_ck(i), h, w, [o_L, o_L], [w/s_L, w/s_L], [], EPD_c{i}, EPD_s);
    
    % Saving results
    res(count,:) = [f_ck(i), h, o_L, s_L, GWP, M_Rd, M_Ed.sm, M_Ed.nc, M_Ed.ns];
    count = count + 1;

end


%% Scenario B
fprintf('\n\n\n---------- Scenario B ---------- \n');

% Printing header
fprintf('\nChange of wall thickness: (C16)\n');

% Characteristic compressive strength of concrete
% (16, 20, and 35 is currently available)
f_ck = 16;   EPD_c = 'c16F';

% Wall thickness
h = [220 180 140 100];

% Diameter of longitunal reinforcement
o_L = 6;

% Spacing between longitunal reinforcement
s_L = 150;

% Looping through each case
for i = 1:size(h, 2)

    % Printing header
    fprintf('\n(C%.f/h%.f/ø%.f/s%.f):   ', f_ck, h(i), o_L, s_L)
    
    % Calculations
    [M_Ed, M_Rd, GWP] = mom_dem(f_ck, h(i), w, [o_L, o_L], [w/s_L, w/s_L], [], EPD_c, EPD_s);
    
    % Saving results
    res(count,:) = [f_ck, h(i), o_L, s_L, GWP, M_Rd, M_Ed.sm, M_Ed.nc, M_Ed.ns];
    count = count + 1;

end


%% Scenario C
fprintf('\n\n\n---------- Scenario C ---------- \n');

% Printing header
fprintf('\nChange of wall thickness with greater precision: (C16)\n');

% Characteristic compressive strength of concrete
% (16, 20, and 35 is currently available)
f_ck = 16;   EPD_c = 'c16F';

% Wall thickness
h = [143 140 137 134 131 128];

% Diameter of longitunal reinforcement
o_L = 6;

% Spacing between longitunal reinforcement
s_L = 150;

% Looping through each case
for i = 1:size(h, 2)

    % Printing header
    fprintf('\n(C%.f/h%.f/ø%.f/s%.f):   ', f_ck, h(i), o_L, s_L)
    
    % Calculations
    [M_Ed, M_Rd, GWP] = mom_dem(f_ck, h(i), w, [o_L, o_L], [w/s_L, w/s_L], [], EPD_c, EPD_s);
    
    % Saving results
    res(count,:) = [f_ck, h(i), o_L, s_L, GWP, M_Rd, M_Ed.sm, M_Ed.nc, M_Ed.ns];
    count = count + 1;

end


%% Scenario D
fprintf('\n\n\n---------- Scenario D ---------- \n');

% Printing header
fprintf('\nChange of longitudinal reinforcement diameter: (C16/h134)\n');

% Characteristic compressive strength of concrete
% (16, 20, and 35 is currently available)
f_ck = 16;   EPD_c = 'c16F';

% Wall thickness
h = 134;

% Diameter of longitunal reinforcement
o_L = [4 6 8 10 12 14];

% Spacing between longitunal reinforcement
s_L = 150;

% Looping through each case
for i = 1:size(o_L, 2)

    % Printing header
    fprintf('\n(C%.f/h%.f/ø%.f/s%.f):   ', f_ck, h, o_L(i), s_L)
    
    % Calculations
    [M_Ed, M_Rd, GWP] = mom_dem(f_ck, h, w, [o_L(i), o_L(i)], [w/s_L, w/s_L], [], EPD_c, EPD_s);
    
    % Saving results
    res(count,:) = [f_ck, h, o_L(i), s_L, GWP, M_Rd, M_Ed.sm, M_Ed.nc, M_Ed.ns];
    count = count + 1;

end


%% Scenario E
fprintf('\n\n\n---------- Scenario E ---------- \n');

% Printing header
fprintf('\nChange of wall thickness: (C16/ø10)\n');

% Characteristic compressive strength of concrete
% (16, 20, and 35 is currently available)
f_ck = 16;   EPD_c = 'c16F';

% Wall thickness
h = [90 93 95 98 101 104];

% Diameter of longitunal reinforcement
o_L = 10;

% Spacing between longitunal reinforcement
s_L = 150;

% Looping through each case
for i = 1:size(h, 2)

    % Printing header
    fprintf('\n(C%.f/h%.f/ø%.f/s%.f):   ', f_ck, h(i), o_L, s_L)
    
    % Calculations
    [M_Ed, M_Rd, GWP] = mom_dem(f_ck, h(i), w, [o_L, o_L], [w/s_L, w/s_L], [], EPD_c, EPD_s);
    
    % Saving results
    res(count,:) = [f_ck, h(i), o_L, s_L, GWP, M_Rd, M_Ed.sm, M_Ed.nc, M_Ed.ns];
    count = count + 1;

end


%% Scenario F
fprintf('\n\n\n---------- Scenario F ---------- \n');

% Printing header
fprintf(['\nChange of spacing betweeen longitudinal reinforcement:' ...
    ' (C16/ø10/h93)\n']);

% Characteristic compressive strength of concrete
% (16, 20, and 35 is currently available)
f_ck = 16;   EPD_c = 'c16F';

% Wall thickness
h = 93;

% Diameter of longitunal reinforcement
o_L = 10;

% Spacing between longitunal reinforcement
s_L = [100 150 200 250];

% Looping through each case
for i = 1:size(s_L, 2)

    % Printing header
    fprintf('\n(C%.f/h%.f/ø%.f/s%.f):   ', f_ck, h, o_L, s_L(i))
    
    % Calculations
    [M_Ed, M_Rd, GWP] = mom_dem(f_ck, h, w, [o_L, o_L], [w/s_L(i), w/s_L(i)], [], EPD_c, EPD_s);
    
    % Saving results
    res(count,:) = [f_ck, h, o_L, s_L(i), GWP, M_Rd, M_Ed.sm, M_Ed.nc, M_Ed.ns];
    count = count + 1;

end


%% Scenario G
fprintf('\n\n\n---------- Scenario G ---------- \n');

% Printing header
fprintf(['\nChange of spacing betweeen longitudinal reinforcement ' ...
    'with greater precision: (C16/ø10/h93)\n']);

% Characteristic compressive strength of concrete
% (16, 20, and 35 is currently available)
f_ck = 16;   EPD_c = 'c16F';

% Wall thickness
h = 93;

% Diameter of longitunal reinforcement
o_L = 10;

% Spacing between longitunal reinforcement
s_L = [110 120 130 140 150 160 170 180];

% Looping through each case
for i = 1:size(s_L, 2)

    % Printing header
    fprintf('\n(C%.f/h%.f/ø%.f/s%.f):   ', f_ck, h, o_L, s_L(i))
    
    % Calculations
    [M_Ed, M_Rd, GWP] = mom_dem(f_ck, h, w, [o_L, o_L], [w/s_L(i), w/s_L(i)], [], EPD_c, EPD_s);
    
    % Saving results
    res(count,:) = [f_ck, h, o_L, s_L(i), GWP, M_Rd, M_Ed.sm, M_Ed.nc, M_Ed.ns];
    count = count + 1;

end


%% Scenario X
fprintf('\n\n\n---------- Scenario X ---------- \n');

% Printing header
fprintf(['\nChange of spacing betweeen longitudinal reinforcement ' ...
    'with greater precision: (C16/ø10/h93)\n']);

% Characteristic compressive strength of concrete
% (16, 20, and 35 is currently available)
f_ck = 16;   EPD_c = 'c16F';

% Wall thickness
h = [90 93 95 98 101 104 107 110 113 116 119];

% Diameter of longitunal reinforcement
o_L = 10;

% Spacing between longitunal reinforcement
s_L = 150;

% Looping through each case
for i = 1:size(h, 2)

    % Printing header
    fprintf('\n(C%.f/h%.f/ø%.f/s%.f):   ', f_ck, h(i), o_L, s_L)
    
    % Calculations
    [M_Ed, M_Rd, GWP] = mom_dem(f_ck, h(i), w, [o_L, o_L], [w/s_L, w/s_L], [], EPD_c, EPD_s);
    
    % Saving results
    res(count,:) = [f_ck, h(i), o_L, s_L, GWP, M_Rd, M_Ed.sm, M_Ed.nc, M_Ed.ns];
    count = count + 1;

end



% Printing conclusion
fprintf(['\n\n\n---------------------------------------- NOTE ', ...
    '---------------------------------------- \n\n'])

% Printing explanation for low design value of bending moment capacity
fprintf(['The estimated design value of the applied bending moment ' ...
    'calculated with the nominal stiffness method, M_Ed_ns,\n' ...
    'decreases greatly in scenario D where the diameter is increased.\n' ...
    'High geometric reinforcement ratio, rho >= 0.01, no additional ' ...
    'strenght can be has from the steel. \n']);


% Specifying column names
colNames = ["f_ck", "h", "o_l", "s_L", "GWP", "M_Rd", "M_Ed_sm", ...
    "M_Ed_nc", "M_Ed_ns"];

% Creates table
tab = table(res(:,1), res(:,2), res(:,3), res(:,4), res(:,5), res(:,6), res(:,7), res(:,8), res(:,9), 'VariableNames', colNames);

% Saving table
writetable(tab, 'wall_analysis.xls')


%--------------------------------------------------------------------------
% Calculation of design values for applied moment and moment capacity

function [M_Ed, M_Rd, GWP] = mom_dem(f_ck, h, w, o_L, n_L, d_L, EPD_c, EPD_s)

    % Definitions
    M = def_w(f_ck, h, w, o_L, n_L, d_L);

    % Geometrical imperfection
    M.e_0 = 5 *10^(-3); %geo_imperfections(0, M);
    
    % First order bending moment
    M = bending_moment(0, M);

    % Effects of creep
    M.varphi_ef = effective_creep(0, M);

    % Moment capacity
    [M_Rd, ~] = moment_capacity(0, M);

    % Applied second order moment
        
        %Simplified method II
        M_Ed.sm = simplified_method(0, M);

        % Nominal curvature
        M_Ed.nc = nominal_curvature(0, M);

        % Nominal stiffness
        [M_Ed.ns, ~] = nominal_stiffness(0, 1.15, M);


    % Global warming potential
    GWP = calc_GWP(0, {EPD_c, EPD_s}, M);

    % Prints results
    fprintf('GWP = %.2f [kg CO2-ep] | ', GWP)
    fprintf('M_Rd = %.2f [kNm] | ', M_Rd*10^(-3))
    fprintf('M_Ed_sm = %.2f [kNm] | ', M_Ed.sm*10^(-3))
    fprintf('M_Ed_nc = %.2f [kNm] | ', M_Ed.nc*10^(-3))
    fprintf('M_Ed_ns = %.2f [kNm] | ', M_Ed.ns*10^(-3))


% Ending function
end