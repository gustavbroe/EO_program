function BL = baseline()

%--------------------------------------------------------------------------
% This function simply returns baseline values for GWP and M_Rd for
% original member SE01 and VES01.
%
% The values are hard-coded on purpose to monitor effects of changes made
% to the core-functions.
%
% Input:    None.
% Output:   A structure array with baseline values.
%--------------------------------------------------------------------------

% Wall, VES01

    % Global Warming Potential (EN 15804 + A1)
    BL.wall.GWP = 1.344821658455886e+03;

    % Moment capacity
    BL.wall.M_Rd = 1.704680049019222e+05;

% Column, SE01

    % Global Warming Potential (EN 15804 + A1)
    BL.column.GWP = 1.768911686694422e+02;

    % Moment capacity
    BL.column.M_Rd = 1.348868769305397e+05;

end