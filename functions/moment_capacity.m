function [M_Rd, rs]  = moment_capacity(info_lev, M)

%--------------------------------------------------------------------------
% This function calculates the design value of the moment capacity, M_Rd 
% for a given member. 
% Accounting for steel's non-linear material stress-strain curve.
%
% Assumptions
%  - The member breaks when the outmost compression fibre is subject to the 
%    ultimate compressive strain of concrete, varepsilon_cu3. (Bi-lienar)
%
% Input:    Information about the member (structure array).
% Output:   Moment capacity in Nm and a number resembling its state.
%--------------------------------------------------------------------------

% Calculcating strength parameters for rectangular stress distribution
P = strength_parameters(M);

% Area of reinforcement steel
P.A_sr = M.n_L .* pi .* (M.o_L/2).^2;                           %[m^2]

% Counter for number of valid calculations
k = 0;                                                          %[-]


% Creating a for-loop that calcaltes the depth of neutral, x, for each
% posible number of rows in compression
for i = 0:size(M.d_L ,2)

    % Index of rows in...

        % compression
        P.index_sc = 1:i;                                       %[-]

        % tension
        P.index_s = i+1:size(M.d_L ,2);                         %[-]

    % First guess at the height of the neutral axis is 100% and 60% of the
    % cross sectional height, h.
    x = [0.80*M.h, 0.70*M.h];                                   %[m]
    
    % Starting counter and setting the change
    % (The initial value of the change has to greater than the limit)
    j = 2;   err = 1;                                           %[-/m]
    
    % Repeating the estimation process until the desired precision is
    % achieved
    while err > 0.0001*10^(-3)
    
        % Change factor (secant methods)
        cf = (x(j-1) - x(j)) /  ... 
            (soe(M, P, x(j-1)).lhs - soe(M, P, x(j)).lhs);      %[-]

        % Improved x-value
        x(j+1) = x(j) - soe(M, P, x(j)).lhs * cf;               %[m]
    
        % Absolute change in approximation of result
        err = abs(x(j) - x(j+1));                               %[m]
    
        % Adding to counter
        j = j + 1;                                              %[-]
        
    end
    
    % Number of steel reinforcement rows with negative strain
    % (ns = negative strain)

        % Compression
        ns_sc = sum(soe(M, P, x(end)).varepsilon_sc < 0);       %[-]

        % Compression
        ns_s = sum(soe(M, P, x(end)).varepsilon_s < 0);         %[-]

    % Storing number of rows in compression for which the calculation is
    % valid
    if (ns_sc + ns_s) == 0
        k = k + 1;
        v_index(k) = i;
        v_x(k) = x(end);
    end

    % Printing infomation about current loop
    if info_lev == 3
        PLI(i, x, soe(M, P, x(end)))
    end

% Ending for-loop
end

% Hopefully there is not two valid solution, but I have to check
if size(v_index, 2) ~= 1
    fprintf('[ERROR] More than one valid solution (%.f)',size(v_index, 2))
    return
end

% Creating index for rows in compression and tension from the valid
% solution found above.
P.index_sc = 1:v_index;   P.index_s = v_index+1:size(M.d_L  ,2);

% Using the same function but using the forces and not the left-hand side.
R = soe(M, P, v_x); 

% Initializing moment capacity (for rectangular stress distribution)
k = 0.5;                                                        %[-]
M_Rd = R.F_c*(M.h/2 - k*P.lambda*v_x);                          %[Nm]

% Added compression row contribution
M_Rd = M_Rd + sum( R.F_sc(P.index_sc) .*  ...
    abs( M.h/2 - M.d_L(P.index_sc) ) );                         %[Nm]

% Added tension row contribution
M_Rd = M_Rd + sum( R.F_s(P.index_s-v_index) .* ... 
    abs( M.h/2 - M.d_L(P.index_s) ) );                          %[Nm]


% Check the reinforcement state of the cross section, by comparing the
% strain in the last row of reinforcement to the ultimate and yeilding.
if R.varepsilon_s(end) > M.varepsilon_uk
    % Meaning that the reinforcement will split/break.
    RS = 'UNDER';   rs = -1;

elseif R.varepsilon_s(end) < M.varepsilon_yd
    % The steel will not yeild/stretch before the concrete breaks
    RS = 'OVER';   rs = 1;

else
    % Other wise the steel is yeilding when the concrete breaks.
    RS = 'NORMAL';   rs = 0;
end



%--------------------------------------------------------------------------
%% Outputting requested information

if info_lev == 1
    fprintf('M_Rd = %.4f [kNm] | ', M_Rd*10^(-3))
    fprintf('x = %.2f [mm] | ', v_x*10^(3))
    fprintf('n_sc = %.f [-] | ', v_index)
    fprintf('%s | ',RS)

elseif info_lev == 2
    fprintf(['The effective depth of the nautral axis is: ' ...
        'x = %.2f [mm]\n'], v_x*10^(3));
    fprintf('The number of rows in compression is: n_sc = %.f\n', v_index);
    fprintf(['The design value of the moment capacity is: ' ...
        'M_Rd = %.2f [kNm] \n'], M_Rd*10^(-3));

    % Prints reinforcement state
    fprintf('The reinforcement state of the cross section is: %s', RS);
    if strcmp(RS, 'UNDER')
        fprintf('\n(varepsilon_s(end) > varepsilon_uk) = (%.2e > %.2e) \n'...
            , R.varepsilon_s(end), M.varepsilon_uk);

    elseif strcmp(RS, 'OVER')
        fprintf('\n(varepsilon_s(end) < varepsilon_yd) = (%.2e < %.2e) \n'...
            , R.varepsilon_s(end), M.varepsilon_yd);

    elseif strcmp(RS, 'NORMAL')
        fprintf(['\n(varepsilon_yd < varepsilon_s(end) < varepsilon_uk) ' ...
            '= (%.2e < %.2e < %.2e) \n'] ...
            , M.varepsilon_yd, R.varepsilon_s(end), M.varepsilon_uk);

    end

end

% Ending function
end



%--------------------------------------------------------------------------
%% State of equilibrium in the direction opposite of the axial load

% The function calculates the force contributions which should equal zero,
% this is done for a given depth of the neutral axis, x.
%
% The equation is the following
% N_Ed = F_c + F_sc - F_s
%
% Rearranging to equla zero
% F_c + F_sc - F_s - N_Ed = 0

function R = soe(M, P, x)

% Strains in the reinforcement rows

    % Compression                                               %[-]
    if P.index_sc ~= 0
        R.varepsilon_sc = P.varepsilon_cu3 .* (x - M.d_L(P.index_sc)) / x;
    else
        R.varepsilon_sc = 0;
    end

    % Tension                                                   %[-]
    if P.index_s ~= 0
        R.varepsilon_s = P.varepsilon_cu3 .* (M.d_L(P.index_s) - x) / x;
    else 
        R.varepsilon_s = 0;
    end


% Stresses in rows
% (The strains can not be greater than the yeild strength of steel)

    % Compression
    R.sigma_sc = min(R.varepsilon_sc * M.E_s, M.f_yd);          %[Pa]

    % Tension
    R.sigma_s = min(R.varepsilon_s * M.E_s, M.f_yd);            %[Pa]


% Resulting force from the cross section components

    % Square stress distribution from concrete
    R.F_c = P.lambda*x*M.w*M.f_cd*P.eta;                        %[N]

    % Compression steel
    R.F_sc = P.A_sr(P.index_sc) .* R.sigma_sc;                  %[N]

    % Tension steel
    R.F_s = P.A_sr(P.index_s) .* R.sigma_s;                     %[N]

% Left hand side of state of equilibrium equation
R.lhs = R.F_c + sum(R.F_sc) - sum(R.F_s) - M.N_Ed;              %[N]

end



%--------------------------------------------------------------------------
%% Additional iformation about the loops

% Print Loop Information
function PLI(i, x, R)

% Header
fprintf('------------- %.f row(s) in compression ------------- \n\n', i);


% Iterations of x
fprintf('Iterations of the effective depth of the neutral axis, x: [mm]\n')
fprintf('%.2f   ', x*10^(3)); fprintf('\n\n');


% Strain
fprintf(['Strain in the rows for the final x-value, varepsilon:' ...
    ' [%%] (from top to bottom)\n']);
if R.varepsilon_sc ~= 0
    fprintf('%.4f   ', R.varepsilon_sc*10^2); 
end

if R.varepsilon_s ~= 0
    fprintf('%.4f   ', R.varepsilon_s*10^2); 
end
fprintf('\n\n');


% Forces
fprintf(['Force from the square effective concrete area: ' ...
    'F_c = %.2f [kN] \n'], R.F_c*10^(-3));
fprintf('Forces from the reinforcement rows, F_s/F_sc: ');
if R.F_sc ~= 0
    fprintf('   %.4f', R.F_sc*10^(-3)); 
end

if R.F_s ~= 0
    fprintf('   %.4f', R.F_s*10^(-3)); 
end
fprintf(' [kN] \n\n\n');


end