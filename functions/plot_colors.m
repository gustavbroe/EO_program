function PC = plot_colors()

%--------------------------------------------------------------------------
% Prints the figures which are currently open, excluding any "Input check"
% windows.
%
% Input:    None
% Output:   Structured array with colers (HEX-code)
%--------------------------------------------------------------------------

% Defining the color of

% Simplified method
    PC.sm = "#FFEC5C";      % Yellow    / 255, 236, 92

    % Nominal curvature
    PC.nc = "#F09041";      % Orange    / 240, 144, 65

    % Nominal stiffness
    PC.ns = "#7CC6FF";      % Blue    / 124, 198, 255
    
    % Moment capacity
    PC.mc = "#D55672";    % Pink      / 213, 86, 114

    % Global warming potential
    PC.gwp = "#1A12FF";     % Purple Blue / 26, 18, 255

    % The color gray
    PC.gray = "#262626";    % Gray      / 38, 38, 38

    % The color red
    PC.red = "#ff6670";     % Red       / 255, 102, 112

    % Yellow
    PC.yellow = "#fff069";  % Yellow    / 255, 240, 105

    % The color green
    PC.green = "#88ff63";   % Green     / 136, 255, 99

end