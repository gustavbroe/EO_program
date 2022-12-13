function plot_cross_section(pf, M, margin)

%--------------------------------------------------------------------------
% Plots the defined cross section in the struc. The function works as an
% control of the input parameters.
%
% The longitunal reinforcement rows are automatically equally spaced within
% the transversal reinforcement.
%
% Input:    A parent object (Only, Axes, UIAxes, etc.)
%           Information about the member (structure array) 
%           A number for the desired distance around the  concrete edge.
% Output:   None (other than the generated plot)
%--------------------------------------------------------------------------


% Color of plot
pcc = '#4060ff';

% Distance from bar's center to concrete's edge
cc = M.h-M.d;

% Number of longitudinal reinforcement rows
rn = size(M.n_L,2);


% CREATING TYPOLOGI MATRIX

% Initializing typologi matrix and a counter
X = zeros(sum(M.n_L), 2); j = 1;

% The outer loop runs equal to the number of rows
for i = 1:rn

    % Distance from the bottom to the current row
    hh = cc + (M.h-2*cc)/(rn-1) * (rn-i);

    % Looping for each bar in the current row
    for f = 1:M.n_L(i)
        
        % Distance from the left side of cross section to bar
        ww = cc + (M.w-2*cc)/(M.n_L(i)-1) * (f-1);

        % Notes the coordinates for the current bar
        X(j, :) = [ww hh];

        % Adding to the counter
        j = j + 1;

    end
end


% PLOTTING

% Gathering information about the current screen
% si = get(0,'screensize');
%     
%     % Height of toolbar
%     th = si(4)*0.06;
% 
%     %Height ployting window header
%     wh = si(4)*0.18;

% Future graph will be added to the same figure
hold(pf, 'on');

% Highlighting plot boundaries
box(pf, 'on');

%Plotting the center of the reinforcement bars
plot(pf, X(:,1), X(:,2), 'r.');

% Adding labels to the centre of the bars
text(pf, X(:,1)+max(M.o_L), X(:,2), num2str((1:sum(M.n_L))'));

% Setting the limit of the axis and making them equal length
rr = margin*min(M.h, M.w);
axis(pf,'equal'); axis(pf, [-rr M.w+rr -rr M.h+rr]);

% Changing the position of the plot
% set(pf, 'Position', [si(3)*0.5 th si(3)*0.5 si(4)-(th+wh)])

% Setting the labels
xlabel(pf, 'x [m]'); ylabel(pf, 'y [m]');


% DRAWING CIRCLES AROUND THE CENTRE DOTS

% Starting counter
z = 1;

% Creating an circle for each row in the typologi matrix
for i = 1:size(M.n_L, 2)

    % Resolution is the middle number (number of points in circle)
    nn = 0:(2*pi)/100:2*pi;

    for j = 1:M.n_L(i)

        % Plotting the individual circels
        plot(pf, X(z,1)+M.o_L(i)/2*cos(nn) ,X(z,2)+M.o_L(i)/2*sin(nn) ...
            ,'Color', pcc, 'Linestyle', '-', 'LineWidth', 1);
        
        % Addding to counter
        z = z + 1;

    end

end


% DRAWING THE MAIN CROSS SECTION GEOMETRY

% Coordinates for the corners
xx = [0 M.w M.w 0 0];
yy = [0 0 M.h M.h 0];

% Plotting a line between the corner points
plot(pf, xx, yy, 'Color', pcc, 'Linestyle', '-', 'LineWidth', 2)


% DRAWING TRANSVERSAL REINFORCEMENT

% Order of point
oop = [1 M.n_L(1) sum(M.n_L) sum(M.n_L)-M.n_L(end)+1 1]';

% Isolating the end bars from the typologi matrix
iX = X(oop, :);

% Drawing the lines
for i = 1:4

    % Offset of the lines in both y- and y-directions for inner and outer
    ios = [0 max(M.o_L) 0 -max(M.o_L) ; max(M.o_L) 0 -max(M.o_L) 0]'/2;
    oos = ios + [0 M.o_T 0 -M.o_T ; M.o_T 0 -M.o_T 0]';
    
    % Points for the first quarter circle
    qcp = pi/2:pi/50:pi;

    % Plotting each INNER lines
    plot(pf, iX(i:i+1 ,1) + ios(i,1), iX(i:i+1 ,2) + ios(i,2), ...
        'Color', pcc, 'Linestyle', '-', 'LineWidth', 1)

    % Plotting each OUTER lines
    plot(pf, iX(i:i+1 ,1) + oos(i,1), iX(i:i+1 ,2) + oos(i,2), ...
        'Color', pcc, 'Linestyle', '-', 'LineWidth', 1)

    % Plotting each inner quarter circle
    plot(pf, max(M.o_L)/2*cos(qcp - pi/2*(i-1)) + iX(i,1), ...
        max(M.o_L)/2*sin(qcp - pi/2*(i-1)) + iX(i,2), ...
        'Color', pcc, 'Linestyle', '-', 'LineWidth', 1)

    % Plotting each inner quarter circle
    plot(pf, (max(M.o_L)/2+M.o_T)*cos(qcp - pi/2*(i-1)) + iX(i,1), ...
        (max(M.o_L)/2+M.o_T)*sin(qcp - pi/2*(i-1)) + iX(i,2), ...
        'Color', pcc, 'Linestyle', '-', 'LineWidth', 1)

end

% Finishing/concluding the plot
hold(pf, 'off') 