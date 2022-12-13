%% Initializing MATLAB
clear; clc; close all hidden;

% Adding function folder to path
addpath('functions')

% Loading climate information and plot colors
CI = material_info(1);   PC = plot_colors();

%--------------------------------------------------------------------------
% User-input (Copied from a_c.m, )

% Array with number of points
np_min = 20;   np_max = 320;
np_points = np_min:10:np_max;

% Characteristic compressive strength of concrete
% (16, 20, 25, 30 [F/R]  |  35, 40 [L/R])
f_ck = 16;   EPD_c = 'c16F';   EPD_s = 'rp';

% Width of cross section
w = 3600;           %[mm]

% Spacing between longitudinal reinforcement
s_L = 150;          %[mm]

% Limits for thickness and longitunal reinforcement diameter
o_L_min = 6;                o_L_max = 16;               %[mm]
h_min = 80;                 h_max = 250;                %[mm]

% Printing of figures
printing = 0;


%--------------------------------------------------------------------------
% Calculating applied moment, moment capacity and GWP for each combination

% Number of longitunal reinforcement bars in each rows
n_L = [w/s_L, w/s_L];

% Defining text for each method
tx = {'Simplified method', 'Nonimal curvature', 'Nominal stiffness'};

% Allocating space for results
res = zeros(size(np_points,2), 3);

% Looping through diffrent number of points on each axes
parfor i = 1:size(np_points, 2)

% Points in intervals
o_L = linspace(o_L_min, o_L_max, np_points(i));
h = linspace(h_min, h_max, np_points(i));

% Rreallocating arrays and starting counter
hh = zeros(np_points(i), np_points(i));   oo=hh;   MM=hh;   rr=hh;
SM = hh;   NC = hh;   NS = hh;   nn = hh;   GG = hh;

% Loop for diameter of longitunal reinforcement
for x = 1:size(o_L, 2)

    % Loop for side length
    for y = 1:size(h, 2)

        % Definitions
        M = def_w(f_ck, h(y), w, [o_L(x), o_L(x)], n_L, []);

        % Geometrical imperfection
        M.e_0 = geo_imperfections(0, M);
        
        % First order bending moment
        M = bending_moment(0, M);
    
        % Effects of creep
        M.varphi_ef = effective_creep(0, M);

        % Moment capacity
        [MM(x,y), rr(x,y)] = moment_capacity(0, M);

        % Applied second order moment
            
            %Simplified method II
            SM(x,y) = simplified_method(0, M);

            % Nominal curvature
            NC(x,y) = nominal_curvature(0, M);

            % Nominal stiffness
            [NS(x,y), nn(x,y)] = nominal_stiffness(0, 1.15, M);

        % Global warming potential
        GG(x,y) = calc_GWP(0, {EPD_c, EPD_s}, M);

        % Storing result
        hh(x,y) = h(y);   oo(x,y) = o_L(x);

    end

end


%--------------------------------------------------------------------------
% Finding the valid cross sections just within the limit of each method


% Plotting the planes and intersection line between the three planes and MM
d = {SM, NC, NS};  vr = {0,0,0};
for k = 1:3

    % Logical array of valid cross sections with adequate moment cap.
    vcs = MM > d{k};

    % Looping through the columns of the logical array
    for f = 1:size(vcs, 2)

        % Checking if there are more than 1 1 in the column
        if ~isempty(find(vcs(:,f) == 1, 1))

            % Index vector for the rows equal 1
            j = find(vcs(:,f) == 1);
    
            % Changen all but the first 1 to 0
            vcs(j(2:end),f) = false;

        end

    % Ending loop
    end

    % Removing the that follow the edge and not the intersection
    j = find(vcs(1,:)==1);   vcs(1,j(2:end)) = false;

    % Storing the results
    vr{k} = vcs;

% Ending for-loop
end

% Storring the results (outside the loop)
res(i, :) = [min(GG(vr{1})), min(GG(vr{2})), min(GG(vr{3}))];

%--------------------------------------------------------------------------
% Outputting results for the current loop

% Printing header
fprintf('\nThe lowest GWP cross section: (%g/%g points) \n', np_points(i), np_points(i)^2)

% Finding index from highest to lowest method
[~, ii] = sort([res(i, :)]);
for z = ii

    % Outputting the lowest possible GWP
    [row, col] = find(GG == min(GG(vr{z})));
    fprintf(' - %s = %.2f [kg CO2-ep] ', tx{z}, GG(row,col));
    fprintf('| s = %.2f [mm] / o_L = %.2f [mm] \n', ...
        hh(row,col), oo(row,col));

end

% Ending the number of points for-loop
end


%--------------------------------------------------------------------------
% Plotting the error as a function of number of points

% Creating figure for the future plot
p1 = figure('Name', 'Absolute error', "PaperSize", [21, 14]);

% Making it the current figure and insuring that future enteries will be
% added to it.
figure(p1);   hold on

% Plot colors
plot_colors = [PC.sm, PC.nc, PC.ns];
%               (Pink)     (Blue)     (Orange)

% Calculating the absolute error
% (Further work: Calculate the ANALYTICAL solution for each method)
abs_error = abs(res - res(end, :));             %[kg CO2-ep]

% The relative error
rel_error = abs_error./res(end, :);             %[-]

% Loop through each method
for i = 1:size(tx,2)

    % Plotting the relative error for each method
    plot(np_points.^2, rel_error(:,i)*100, "Color", plot_colors(i), ...
        "DisplayName", tx{i}, "Marker", "o", "MarkerSize", 9, ...
        "LineStyle", "-", "LineWidth", 2);

end

% Changing plot legends
xlabel('Number of points, {\it n_{points}} [-]')
ylabel('Relative error, {\it error_{rel}} [%]')

% Concluding creation of the plot
conclude_plot(p1, 'right');


%--------------------------------------------------------------------------
% Plotting the error as a function of number of points

% Creating figure for the future plot
p2 = figure('Name', 'Error', "PaperSize", [21, 14]);

% Making it the current figure and insuring that future enteries will be
% added to it.
figure(p2);   hold on

% Calculating the error
% (Further work: Calculate the ANALYTICAL solution for each method)
error = res - res(end, :);             %[kg CO2-ep]

% Loop through each method
for i = 1:size(tx,2)

    % Plotting the relative error for each method
    plot(np_points.^2, error(:,i), "Color", plot_colors(i), ...
        "DisplayName", tx{i}, "Marker", "o", "MarkerSize", 9, ...
        "LineStyle", "-", "LineWidth", 2);

end

% Changing plot legends
xlabel('Number of points, {\it n_{points}} [-]')
ylabel('Error, {\it error} [%]')

% Concluding creation of the plot
conclude_plot(p2, 'left');


%--------------------------------------------------------------------------
% Printing

if printing == 1

% Path to the results folder
path_folder = '.\results\full_analysis_column';

% Number of figures
number_figures = 1;

% Name of the figures
name_file = append('w_np', string(np_min), '_', string(np_max), ...
    '_nL', join(string([2,2]),'_'), '_C', sprintf('%.f_',f_ck), ...
    'Convergence');

% Printing both a PDF and EPS
print_figures(path_folder, name_file, number_figures);

end


%--------------------------------------------------------------------------
% Concluding plot

function conclude_plot(plot_identifier, screen_side)

% Adding grid lines
grid on

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

% Concluding creation of the plot
hold off

end
