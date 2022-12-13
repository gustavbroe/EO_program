%% Initializing Matlab
clear, clc, close;

% Adding function folder to path
addpath('functions')

% Starts copy of command line output to .txt file
start_diary("c_demonstrate_Yasemin.txt")


%% Scenario 1
fprintf('\n---------- Scenario 1 ---------- \n');

% Printing header
fprintf('Change the concrete strength class: \n')

% Stregnths classes
f_ck = [25 30 35 40];

% Loop for each stregnth class
for i = f_ck

    % Printing header
    fprintf('\nC%.f:   ', i)

    % Calculations
    [M_Ed, M_Rd] = mom_dem(i, 240, 360, [16 16 16], [5 2 5], [150]);

end


% Printing header
fprintf('\n\nChange the Diameter of longitudinal reinforcement (C25): \n')

% Diameter of longitudinal reinforcement
oL = zeros(4,3) + [16 20 25 32]';

% Loop for each diameter of longitudinal reinforcement
for i = 1:size(oL, 1)

    % Printing header
    fprintf('\nØ%.f:   ', oL(i,1))

    % Calculations
    [M_Ed, M_Rd] = mom_dem(25, 240, 360, oL(i,:), [5 2 5], [150]);

end



%% Scenario 2
fprintf('\n\n\n---------- Scenario 2 ---------- \n');

% Printing header
fprintf(['\n\nChange the height, width, and diameter of longitudinal ' ...
    'reinforcement (C25): \n']);

% Diameter of longitudinal reinforcement
oL = zeros(4,2) + [16 20 25 32]';

% Height and width of the cross sections
h = [305 280 250 225];
w = h;

% Loop for each diameter of longitudinal reinforcement
for i = 1:size(oL, 1)

    % Printing header
    fprintf('\n(h, w, ø) = (%.f, %.f, %.f):   ', h(i), w(i), oL(i,1))

    % Calculations
    [M_Ed, M_Rd] = mom_dem(25, h(i), w(i), oL(i,:), [2 2], []);

end



%% Scenario 3
fprintf('\n\n\n---------- Scenario 3 ---------- \n');

% Printing header
fprintf(['\n\nChange the height, width, and diameter of longitudinal ' ...
    'reinforcement (C25): \n']);

% Diameter of longitudinal reinforcement
oL = zeros(4,2) + [16 20 25 32]';

% Height and width of the cross sections
h = [280 250 225 215];
w = h;

% Loop for each diameter of longitudinal reinforcement
for i = 1:size(oL, 1)

    % Printing header
    fprintf('\n(h, w, ø) = (%.f, %.f, %.f):   ', h(i), w(i), oL(i,1))

    % Calculations
    [M_Ed, M_Rd] = mom_dem(25, h(i), w(i), oL(i,:), [3 3], []);

end


%% Scenario 4
fprintf('\n\n\n---------- Scenario 4 ---------- \n');

% Printing header
fprintf(['\n\nChange the height, width, and diameter of longitudinal ' ...
    'reinforcement (C25): \n']);

% Diameter of longitudinal reinforcement
oL = zeros(4,3) + [16 20 25 32]';

% Height and width of the cross sections
h = [275 250 230 215];
w = h;

% Loop for each diameter of longitudinal reinforcement
for i = 1:size(oL, 1)

    % Printing header
    fprintf('\n(h, w, ø) = (%.f, %.f, %.f):   ', h(i), w(i), oL(i,1))

    % Calculations
    % (For this to be correct the diameter of the transversal
    % reinforcement has to be neglected.)
    [M_Ed, M_Rd] = mom_dem(25, h(i), w(i), oL(i,:), [3 2 3], [h(i)-75-15]);

end

% Printing conclusion
fprintf(['\n\n\n---------------------------------------- NOTE ', ...
    '---------------------------------------- \n\n'])

% Printing explanation for low design value of bending moment capacity
fprintf(['The estimated applied moment from the simplified method is ' ...
    'around two times greater than the actual applied moment.\n' ...
    'This is caused by the assumption of the column having a free end' ...
    ' and therfore not being simply supported, which doubles the ' ...
    'effect length of the column, \n' ...
    'resulting in a quadrupling of the estimated deflection. \n']);


% Ends copy of command line output to .txt file
diary off



%--------------------------------------------------------------------------
% Calculation of design values for applied moment and moment capacity

function [M_Ed, M_Rd] = mom_dem(f_ck, h, w, o_L, n_L, d_L)

    % Loading definitions 
    M = def_Yasemin(f_ck, h, w, o_L, n_L, d_L);
    
    % Geometrical imperfection
    M.e_0 = 5*10^(-3); %geo_imperfections(M);
    
    % First order bending moment
    M = bending_moment(0, M);

    % Effects of creep
    M.varphi_ef = effective_creep(0, M);
    
    % Checking if second order effects can be neglected
    soe = check_soe(0, M);

    % Method of nominal stiffness
    M_Ed.ns = nominal_stiffness(0, 1.5, M);

    % Method of nominal curvature
    M_Ed.nc = nominal_curvature(0, M);

    % Simplified method II
    M_Ed.sm = simplified_method(1, M);

    % Moment capacity
    M_Rd = moment_capacity(1, M);


% Ending function
end


%--------------------------------------------------------------------------
% Enable new diary 

function start_diary(diary_name)

% Check if a previous diary exists
if isfile(diary_name)

    % Ensuring that deleted files are moved to the Recycling bin
    recycle('on');

    % Deletes the previous diary
    delete(diary_name)

end

% Turns on diary
diary (diary_name)


% Ending function
end


%--------------------------------------------------------------------------
% Nested definition function

function M = def_Yasemin(f_ck, h, w, o_L, n_L, d_L)

%--------------------------------------------------------------------------
% This fanction is defined be the user and converts the input to a
% structual array carrying all the information.
%
% Input:    Height and width of the cross section.
% Output:   A structure array with geometrical and mechanical properties.
%--------------------------------------------------------------------------


%% Input parameters [SI-Units]
% ASSUMPTIONS
    % Concrete age at loading
    % (In-situ concrete has to remain in the cast for at less 3 days.
    % Lisbet Sand, DTU Campusservice, said during the tour of B112, that
    % the each level was constructed over 1-2 weeks.
    % In-situ concrete has to remain in the cast for at less 3 days,
    % then the rest of building is built in the following 2-4 month.
    % 30 days is assumed, as load is progressively added.)
    M.t0 = 30;                         %[days]

    % Eccentricity of load
    % (A good assumption, advised by PG)
    M.e_1 = 55 *10^(-3);                        %[mm]

    % The characterustuc strain at maximum force
    % (For Class B reinforcement steel, smalles value selected)
    M.varepsilon_uk = 5.0 /100;

% CROSS SECTION
    % Height and width
    M.h = h *10^(-3);                           %[mm]
    % Width
    M.w = w *10^(-3);                           %[mm]
    % Area
    M.A_c = M.h * M.w;                          %[m^2]
    % Moment of inertia (rectangular)
    M.I_c = 1/12*M.w*M.h^3;                     %[m^4]

% REINFORCEMENT
    % Thickness of concret cover
    M.c = 15 *10^(-3);                          %[mm]
    % Diameter (longitudinal)
    M.o_L = o_L *10^(-3);                   %[mm]
    % Number of bars in each row [no less than 2]
    M.n_L = n_L;                              %[-]
    % Diameter (transversal)
    M.o_T = 6 *10^(-3);                         %[mm]
    % Depth to first row of longitudinal bars
    M.d_sc = M.c + M.o_T + M.o_L(1)/2;        %[m]
    % Effective depth
    M.d = M.h - (M.c + M.o_T + M.o_L(end)/2);    %[m]
    % Depth of the bars center of each row
    M.d_L = [M.d_sc, (d_L*10^(-3)), M.d];      %[m]
    % Number of reinforcement rows in compression (First estimation)
    M.n_sc = 1;

% COLUMN
    % Length
    M.L = 4.180;                                 %[m]
    % Effective length factor
    M.L_0 = 2*M.L;                              %[m]

% PARTIAL FACTORS
    % Control class (Normal)
    M.gamma_3 = 1.00;                           %[-]

% CONCRETE
    % Characteristic compressive strength (C30)
    M.f_ck = f_ck *10^6;                          %[MPa]
    % Secant modulus at 40% of the average compressive strength
    M.E_cm = 33 *10^9;                          %[GPa]
    % Ultimate strain
    M.varepsilon_cu = 0.35 /100;                %[%]
    % Partial factor (Reinforced and under compression)
    M.gamma_c = 1.45 * M.gamma_3;               %[-]

% STEEL
    % Modolus of elesticity
    M.E_s = 200 *10^9;                          %[GPa]
    % Yeild strength
    M.f_yk = 550 *10^6;                         %[MPa]
    % Partial factor (Slap armering)
    M.gamma_s = 1.20 * M.gamma_3;               %[-]

% LOAD
    % Design value of axial load (ULS)
    M.N_Ed = 455 *10^3;                         %[kN]
    % Quasi-permanent load (SLS)
    M.N_0qp = 1 *10^3;                        %[kN]
    % Transverse line load (from wind)
    M.p_tv = 0 *10^3;                          %[kN/m]


%% Preliminary calculations

% AREAS
    % Total area of steel
    M.A_s = 0;
    for i = 1:size(M.o_L, 2)
        M.A_s = M.A_s + M.n_L(i) * pi * (M.o_L(i)/2)^2;         %[m^2]
    end
    % Second moment of area (about center)
    M.I_s = 0;
    for i = 1:size(M.o_L, 2)
        % Check if the current row is over or under the centre
        if M.d_L(i) <= M.h/2
            M.I_s = M.I_s + M.n_L(i) * pi * (M.o_L(i)/2)^2 ...
            * (M.h/2 - M.d_L(i))^2;                             %[m^2]
        else
            M.I_s = M.I_s + M.n_L(i) * pi * (M.o_L(i)/2)^2 ...
            * (M.d_L(i) - M.h/2)^2;                             %[m^2]
        end
    end

% CONVERTING CHARACTERISTIC VALUES TO DESIGN VALUES
    % Yeild strength of steel
    M.f_yd = M.f_yk / M.gamma_s;                %[Pa]
    % Compressive strength of concrete
    M.f_cd = M.f_ck / M.gamma_c;                %[Pa]
    % Modulus of concrete
    M.E_cd = M.E_cm / M.gamma_c;                %[Pa]
    % Yield strain of steel
    M.varepsilon_yd = M.f_yd/M.E_s;             %[-]
    % strain at maximum force of steel
    M.varepsilon_ud = M.varepsilon_uk/M.gamma_s;%[-]

% REINFORCEMENT RATIO
    % Geometric reinforcement ratio
    M.rho = M.A_s / M.A_c;                      %[-]
    % Mechanical reinforcement ratio
    M.omega = (M.A_s*M.f_yd) / (M.A_c*M.f_cd);  %[-]

% ADDITIONAL GEOMETRIC PARAMETERS
    % Radius of gyration of the uncracked concrete section
    M.i_c = sqrt(M.I_c / M.A_c);                %[m]
    % Slenderness ratio
    M.lambda = M.L_0 / M.i_c;                   %[-]

% ADDITIONAL GEOMETRIC PARAMETERS
    % Calculating the relative axial force
    M.n = M.N_Ed / (M.A_c*M.f_cd);              %[-]        %{5.8.7.2 (2)}
    
end