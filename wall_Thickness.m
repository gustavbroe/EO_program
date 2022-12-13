%% Initializing MATLAB
clear; clc; close all hidden;

% Adding function folder to path
addpath('functions')

% Loading climate information, plot colors, baseline results
CI = material_info(1);   PC = plot_colors();   BL = baseline_line();

%--------------------------------------------------------------------------
%% User-input

% Characteristic compressive strength of concrete
% (16, 20, 25, 30 [F/R]  |  35, 40 [L/R])
f_ck = 16;   EPD_c = 'c16F';   EPD_s = 'rp';

% Width of wall
w = 3600;                       %[mm]

% The start, interval and end of side lengths
h_min = 100;   h_max = 250;
h = h_min:0.1:h_max;            %[mm]  

% Diameter of longitunal reinforcement
o_L = [6, 6];                   %[mm]

% Number of bars in each row
s_L = 150;                      %[mm]

% Printing
printing = 0;

% Plot GWP on a second axis
PG = false;

% Print info til command line
PTC = 0;

% Input control (enable if desired)
% input_control(0.03, def_w(f_ck, mean(h), w, o_L, [w/s_L, w/s_L], []));

%--------------------------------------------------------------------------
%% Analysis of bending moment as a function of side length

% Preallocating results array
q = 8;   res = zeros(size(h, 2), q);

% The start, interval and end of side lengths
for i = 1:size(h,2)
    
    % Displaying current side lenght
    if PTC == 1
        fprintf('\ns = %.4f | ', h(i))
    end

    % Loading definitions 
    M = def_w(f_ck, h(i), w, o_L, [w/s_L, w/s_L], []);
    
    % Geometrical imperfection
    M.e_0 = geo_imperfections(PTC, M);
    
    % First order bending moment
    M = bending_moment(PTC, M);

    % Effects of creep
    M.varphi_ef = effective_creep(PTC, M);
    
    % Checking if second order effects can be neglected
    soe = check_soe(PTC, M);

    % Simplified method II
    M_Ed_sm = simplified_method(PTC, M);

    % Method of nominal curvature
    M_Ed_nc = nominal_curvature(PTC, M);

    % Method of nominal stiffness
    [M_Ed_ns, nr] = nominal_stiffness(PTC*2, 1.2, M);

    % Moment capacity
    [M_Rd, rs] = moment_capacity(PTC, M);

    % Global warming potential
    GWP = calc_GWP(0, {EPD_c, EPD_s}, M);

    % Saving the results of the current loop
    res(i,:) = [M.h soe M_Ed_sm M_Ed_nc M_Ed_ns, M_Rd, GWP, nr];

end

%% Plotting the results regarding moment
 
% Makersize
MS = 15;

% Creating figure for the future plot
p1 = figure('Name', 'Moment and GWP', "PaperSize", [21, 14]);

% Making it the current figure and insuring that future enteries will be
% added to it.
figure(p1);   hold on

% The results are only plotted if there is more than 50 results
if size(res,1) < 50
    return
end

% Define margins
margin = 0.05;   

% Specifying axis color
if PG == true
    colororder(string({PC.gray, PC.gwp}));
    yyaxis left;
end

% Ensuring further plots are in same figure
hold on;

% Plotting the desgin value of the bending moment from the three methods
    % First order
    plot([min(res(:,1)) ; max(res(:,1))], [M.M_0Ed ; M.M_0Ed], ...
        "Marker", "none", "Color", PC.gray, "LineStyle", "-", ...
        "DisplayName", 'First order moment, {\it M_0_E_d}');
    % Nominal stiffness
    plot(res(:,1), res(:,3), "MarkerSize", MS, "Marker", ".", ...
        "Color", PC.sm, "LineStyle", "none", 'MarkerSize', 10, ...
        "DisplayName", 'Simplified method, {\it M_E_d_s_m}');
    % Nominal curvature
    plot(res(:,1), res(:,4), "MarkerSize", MS, "Marker", ".", ...
        "Color", PC.nc, "LineStyle", "none", 'MarkerSize', 10, ...
        "DisplayName", 'Nominal curvature, {\it M_E_d_n_c}');
    % Simplified nominal curvature method
    plot(res(:,1), res(:,5), "MarkerSize", MS, "Marker", ".", ...
        "Color", PC.ns, "LineStyle", "none", 'MarkerSize', 10, ...
        "DisplayName", 'Nominal stiffness, {\it M_E_d_n_s}');

% Plotting the design value of the moment capacity with and without sc
plot(res(:,1), res(:,6), "MarkerSize", MS, "Marker", ".", ...
    'Color', PC.mc, 'LineStyle', 'none', "DisplayName", ...
    'Moment capacity, {\it M_R_d}', 'MarkerSize', 10);


% Setting the labels
xlabel('Thickness of wall, {\it h} [m]'); 
ylabel('Design value of the moment parameters, {\it M} [Nm]');

% Setting axis
    % Maximum y-value
    ymax = max(max(res(:,3:end-2)));
    % Minimum y-value
    ymin = min(min(res(:,3:end-2)));
    % Height of plot
    hop = max(max(res(:,3:end-2))) - min(min(res(:,3:end-2)));
    % Changing the axis
    axis([res(1,1), res(end,1), ymin - hop*margin, ...
    ymax + hop*margin]);


% Checking if second order effects can be ignored for any side length
if sum(res(:,2)==0) > 0

    % PLotting a vertical line for when the second order effects 
    % can be ignored.
        % Side length at limit
        slim = min(res(res(:,2)==0, 1));
        % Index of slim
        slim_i = find(res(:,1)==slim);
        %Plotting
        plot([slim ; slim], [ymin - hop*margin ; ymax + hop*margin], ...
        "Marker", "none", "Color", PC.gray, "LineStyle", "--", ...
        "HandleVisibility", "off");

end


% Checking if the cross section is non-reinforced for any side length
if sum(res(:,end)==1) > 0

    % PLotting a vertical line for when the second order effects 
    % can be ignored.
        % Side length at limit
        slim = min(res(res(:,end)==1, 1));
        % Index of slim
        slim_i = find(res(:,1)==slim);
        %Plotting
        plot([slim ; slim], [ymin - hop*margin ; ymax + hop*margin], ...
        "Marker", "none", "Color", PC.red, "LineStyle", "-.", ...
        "HandleVisibility", "off");

end

% Adding plot title
title(sprintf(['(N_E_d = %.1f [kN] / f_c_k = %.f [MPa] / L_0 = %.2f [m]' ...
    ' / c+o_T = %.f [mm] / o_L = %.f [mm])'] ...
    , M.N_Ed*10^(-3), M.f_ck*10^(-6), M.L_0, (M.c+M.o_T)*10^(3) ...
    , o_L(1)), "FontSize", 11, "Color", [0.15, 0.15, 0.15], ...
    "FontWeight", "normal" ,"Units", "normalized", ...
    "Position", [0.5, 1.050, 0.5]);

% Concluding plot
conclude_plot(p1, 'left');


%--------------------------------------------------------------------------
%% Plotting the results regarding GWP

if PG == true

% Plotting to the right axis now
yyaxis right

% Plotting teh GWP data
plot(h*10^(-3), res(:,7), "MarkerSize", MS, "Marker", ".", ...
    "Color", PC.gwp, "LineStyle", "none", "DisplayName", ...
    'Global warming potential, {\it GWP}');

% Setting the ylabel
ylabel('Global Warming Potential,{\it GWP} [kg CO_2-eq]');

% Ending the plot
hold off

% Adding legends
legend("Location","northwest");

end


%--------------------------------------------------------------------------
% Printing

if printing == 1

    % Path to the results folder
    path_folder = '.\results\Wall\Two_dimensional';
    
    % Number of figures
    number_figures = 1;
    
    % Name of the figures
    name_file = append('h', string(h_min), '_', string(h_max), ...
        '_nL', join(string(M.n_L),'_'), '_oL', ...
        sprintf('%.f',M.o_L(end)*10^(3)), ...
        '_C', sprintf('%.fh',M.f_ck*10^(-6)));
    
    % Printing both a PDF and EPS
    print_figures(path_folder, name_file, number_figures);

end


%--------------------------------------------------------------------------
% Outputting relevant data

% Side for the slenderness is equal to the limiting value
j = find(res(:,2)==1);
fprintf('s_lim = %.2f [mm] \n\n', res(max(j),1)*1E3);

% Side for when the cross section is non-reinforced
j = find(res(:,end)==1);
fprintf('s_non-reinforced = %.2f [mm] \n\n', res(min(j),1)*1E3);

% Side for the side length for discontinuity of NS
j = find(M.A_s./((w*h)*1E-3).^2 > 0.01);
fprintf('s_ns = %.2f [mm] \n\n', res(max(j),1)*1E3);


%--------------------------------------------------------------------------
% Outputting GWP-GHG values

% Allocating for speed
PR = zeros(4,3);

% Each method
for i = 3:5

    % Finding the point when the capacity is greater then the applied
    % moment
    j = find(res(:,6) > res(:,i), 1, 'first');

    % Creating a column with all the data
    PR(:,i-2)  = [res(j,1)*1E3, res(j,i)*1E-3, res(j,7), ...
        (res(j,7)-BL.wall.GWP)/(BL.wall.GWP)]';

end

% Defining row and column names
RowNames = {'s[mm]', 'M_Rd[kNm]', 'GWP-GHG[kg CO2-eq.]', 'GWP_change[-]'};
ColNames = {'Simplified method', 'Nominal curvature', 'Nominal stiffness'};

% Creating table
temp_tab = table(PR, 'RowNames', RowNames);
res_tab = splitvars(temp_tab);
res_tab.Properties.VariableNames = ColNames;

% Printing table
disp(res_tab);


%--------------------------------------------------------------------------
% Concluding plot

function conclude_plot(plot_identifier, screen_side)

% Adding grid lines
grid on;

% Gathering display information
    % Current screen size
    si = get(0,'screensize');
    % Height of toolbar
    th = si(4)*0.06;

% The height/width-ratio of the figure PaperSize
hwr = plot_identifier.PaperSize(2)/plot_identifier.PaperSize(1);

% Changing the position of the plot
if strcmp(screen_side, 'left')

    % To the left side
    set(plot_identifier, 'Position', ...
        [0 th si(3)*0.5 si(3)*0.5*hwr])

elseif strcmp(screen_side, 'right')

    % To the right side
    set(plot_identifier, 'Position', ...
        [si(3)*0.5 th si(3)*0.5 si(3)*0.5*hwr])

end

% Inserting legends
legend;

% Finding the 'Axes'-objects
ax = findall(plot_identifier, 'Type', 'Axes');

%     % Changing the font name
%     set(ax, 'FontName', 'Segoe UI');

    % Adding minor grid lines for both axis
    set(ax, 'XMinorGrid', 'on');
    set(ax, 'YMinorGrid', 'on');


% Concluding creation of the plot
hold off

end