%% Initializing Matlab
clear, clc, close all hidden;

% Adding function folder to path
addpath('functions')

% Starts copy of command line output to .txt file
start_diary("c_demonstrate_Emilie.txt")

% Width of wall
w = 3600;           %[mm]

%% Scenario 1
fprintf('\n---------- Scenario A ---------- \n');

% Printing header
fprintf('\nChange the concrete strength class: \n')

% Stregnths classes
f_ck = [40 35 30 25 20];

% Wall thickness
h = 220;

% Diameter of longitunal reinforcement
o_L = 6;

% Spacing between longitunal reinforcement
s_L = 150;

% Loop for each scenario
for i = 1:size(f_ck, 2)

    % Printing header
    fprintf('\nC%.f/h%.f/ø%.f/s%.f:   ', f_ck(i), h, o_L, s_L)

    % Calculations
    [M_Ed, M_Rd] = mom_dem(f_ck(i), h, w, o_L, w/s_L*4, []);

end


%% Scenario 2
fprintf('\n\n\n---------- Scenario B ---------- \n');

% Printing header
fprintf('\nChange of concrete strength class and thickness: \n');

% Stregnths classes
f_ck = [40 35 30 25 20];

% Wall thickness
h = [97 98 99 101 104];

% Diameter of longitunal reinforcement
o_L = 6;

% Spacing between longitunal reinforcement
s_L = 150;

% Loop for each scenario
for i = 1:size(f_ck, 2)

    % Printing header
    fprintf('\nC%.f/h%.f/ø%.f/s%.f:   ', f_ck(i), h(i), o_L, s_L)

    % Calculations
    [M_Ed, M_Rd] = mom_dem(f_ck(i), h(i), w, o_L, w/s_L*4, []);

end



%% Scenario 3
fprintf('\n\n\n---------- Scenario C ---------- \n');

% Printing header
fprintf(['\nChanging thickness and spacing between longitunal ' ...
    'reinforcement (Ø5): \n']);

% Stregnths classes
f_ck = 25;

% Wall thickness
h = [100 110 119 126];

% Diameter of longitunal reinforcement
o_L = 5;

% Spacing between longitunal reinforcement
s_L = [100 150 200 250];

% Loop for each scenario
for i = 1:size(h, 2)

    % Printing header
    fprintf('\nC%.f/h%.f/ø%.f/s%.f:   ', f_ck, h(i), o_L, s_L(i))

    % Calculations
    [M_Ed, M_Rd] = mom_dem(f_ck, h(i), w, o_L, w/s_L(i)*4, []);

end


% Intermediate header
fprintf(['\n\n\nChanging thickness and spacing between longitunal ' ...
    'reinforcement (Ø6): \n'])

% Wall thickness
h = [93 101 109 115];

% Diameter of longitunal reinforcement
o_L = 6;

% Loop for each scenario
for i = 1:size(h, 2)

    % Printing header
    fprintf('\nC%.f/h%.f/ø%.f/s%.f:   ', f_ck, h(i), o_L, s_L(i))

    % Calculations
    [M_Ed, M_Rd] = mom_dem(f_ck, h(i), w, o_L, w/s_L(i)*4, []);

end


% Intermediate header
fprintf(['\n\n\nChanging thickness and spacing between longitunal ' ...
    'reinforcement (Ø8): \n'])

% Wall thickness
h = [88 91 96 100];

% Diameter of longitunal reinforcement
o_L = 8;

% Loop for each scenario
for i = 1:size(h, 2)

    % Printing header
    fprintf('\nC%.f/h%.f/ø%.f/s%.f:   ', f_ck, h(i), o_L, s_L(i))

    % Calculations
    [M_Ed, M_Rd] = mom_dem(f_ck, h(i), w, o_L, w/s_L(i)*4, []);

end


% Intermediate header
fprintf(['\n\n\nChanging thickness and spacing between longitunal ' ...
    'reinforcement (Ø10): \n'])

% Wall thickness
h = [93 89 90 93];

% Diameter of longitunal reinforcement
o_L = 10;

% Loop for each scenario
for i = 1:size(h, 2)

    % Printing header
    fprintf('\nC%.f/h%.f/ø%.f/s%.f:   ', f_ck, h(i), o_L, s_L(i))

    % Calculations
    [M_Ed, M_Rd] = mom_dem(f_ck, h(i), w, o_L, w/s_L(i)*4, []);

end


fprintf(['\n\n\n------------------------- ' ...
    'C25 -> C20 ' ...
    '-------------------------\n\n'])

% Printing header
fprintf(['\nChanging thickness and spacing between longitunal ' ...
    'reinforcement (Ø5): \n']);

% Stregnths classes
f_ck = 20;

% Wall thickness
h = [102 112 121 128];

% Diameter of longitunal reinforcement
o_L = 5;

% Spacing between longitunal reinforcement
s_L = [100 150 200 250];

% Loop for each scenario
for i = 1:size(h, 2)

    % Printing header
    fprintf('\nC%.f/h%.f/ø%.f/s%.f:   ', f_ck, h(i), o_L, s_L(i))

    % Calculations
    [M_Ed, M_Rd] = mom_dem(f_ck, h(i), w, o_L, w/s_L(i)*4, []);

end


% Intermediate header
fprintf(['\n\n\nChanging thickness and spacing between longitunal ' ...
    'reinforcement (Ø6): \n'])

% Wall thickness
h = [96 104 111 117];

% Diameter of longitunal reinforcement
o_L = 6;

% Loop for each scenario
for i = 1:size(h, 2)

    % Printing header
    fprintf('\nC%.f/h%.f/ø%.f/s%.f:   ', f_ck, h(i), o_L, s_L(i))

    % Calculations
    [M_Ed, M_Rd] = mom_dem(f_ck, h(i), w, o_L, w/s_L(i)*4, []);

end


% Intermediate header
fprintf(['\n\n\nChanging thickness and spacing between longitunal ' ...
    'reinforcement (Ø8): \n'])

% Wall thickness
h = [94 95 99 103];

% Diameter of longitunal reinforcement
o_L = 8;

% Loop for each scenario
for i = 1:size(h, 2)

    % Printing header
    fprintf('\nC%.f/h%.f/ø%.f/s%.f:   ', f_ck, h(i), o_L, s_L(i))

    % Calculations
    [M_Ed, M_Rd] = mom_dem(f_ck, h(i), w, o_L, w/s_L(i)*4, []);

end


% Intermediate header
fprintf(['\n\n\nChanging thickness and spacing between longitunal ' ...
    'reinforcement (Ø10): \n'])

% Wall thickness (Assuming that the first thickness is 93 [mm])
h = [93 95 96 97];

% Diameter of longitunal reinforcement
o_L = 10;

% Loop for each scenario
for i = 1:size(h, 2)

    % Printing header
    fprintf('\nC%.f/h%.f/ø%.f/s%.f:   ', f_ck, h(i), o_L, s_L(i))

    % Calculations
    [M_Ed, M_Rd] = mom_dem(f_ck, h(i), w, o_L, w/s_L(i)*4, []);

end

% Printing conclusion
fprintf(['\n\n\n---------------------------------------- NOTE ', ...
    '---------------------------------------- \n\n'])

% Printing explanation for low design value of bending moment capacity
fprintf(['The design value of bending moment capacity is ' ...
    'lower than the report for cases which are OVER-reinforced, \n' ...
    'meaning that the strain in the bottom row of reinforcement, ' ...
    'varepsilon_s is greater than the steel''s yeild strain, ' ...
    'varepsilon_yd. \nThe result calculated above assumes that ' ...
    'no additional strength from the steel can be had when the strain' ...
    'exceed the yeild strain.\n']);


% Ends copy of command line output to .txt file
diary off



%--------------------------------------------------------------------------
% Calculation of design values for applied moment and moment capacity

function [M_Ed, M_Rd] = mom_dem(f_ck, h, w, o_L, n_L, d_L)

        % Definitions
        M = def_Emilie(f_ck, h, w, o_L, n_L, d_L);

        % Geometrical imperfection
        M.e_0 = 5 *10^(-3); %geo_imperfections(0, M);
        
        % First order bending moment
        M = bending_moment(1, M);
    
        % Effects of creep
        M.varphi_ef = effective_creep(0, M);

        % Moment capacity
        [M_Rd, ~] = moment_capacity(1, M);

        % Applied second order moment
            
            %Simplified method II
            M_Ed.sm = simplified_method(1, M);

            % Nominal curvature
            M_Ed.nc = nominal_curvature(0, M);

            % Nominal stiffness
            [M_Ed.ns, ~] = nominal_stiffness(0, 1.4, M);


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

function M = def_Emilie(f_ck, h, w, o_L, n_L, d_L)

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
    M.e_1 = 25 *10^(-3);                        %[mm]

    % The characterustuc strain at maximum force
    % (For Class B reinforcement steel, smalles value selected)
    M.varepsilon_uk = 5.0 /100;

    % Concrete age at loading
    % (In-situ concrete has to remain in the cast for at less 3 days.
    % Lisbet Sand, DTU Campusservice, said during the tour of B112, that
    % the each level was constructed over 1-2 weeks.
    % In-situ concrete has to remain in the cast for at less 3 days,
    % then the rest of building is built in the following 2-4 month.
    % 30 days is assumed, as load is progressively added.)
    M.t0 = 30;                         %[days]

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
    M.o_T = 5 *10^(-3);                         %[mm]
    % Depth to first row of longitudinal bars
    if size(M.n_L, 2) > 1
        M.d_sc = M.c + M.o_T + M.o_L(1)/2;      %[m]
    else
        M.d_sc = [];
    end
    % Effective depth
    M.d = M.h - (M.c + M.o_T + M.o_L(end)/2);    %[m]
    % Depth of the bars center of each row
    M.d_L = [M.d_sc, (d_L*10^(-3)), M.d];      %[m]
    % Number of reinforcement rows in compression (First estimation)
    M.n_sc = 1;

% COLUMN
    % Length
    M.L = 4.430;                                 %[m]
    %Effective length factor
    M.L_0 = 1*M.L;                              %[m]

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
    M.E_s = 220 *10^9;                          %[GPa]
    % Yeild strength
    M.f_yk = 550 *10^6;                         %[MPa]
    % Partial factor (Slap armering)
    M.gamma_s = 1.20 * M.gamma_3;               %[-]

% LOAD
    % Design value of axial load (ULS)
    M.N_Ed = 540.3 *10^3;                       %[kN]
    % Quasi-permanent load (SLS)
    M.N_0qp = 1 *10^3;                          %[kN]
    % Transverse line load (from wind)
    M.p_tv = 1.12 *10^3;                          %[kN/m]


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