function soe = check_soe(info_lev, M)

%--------------------------------------------------------------------------
% In accordance with EC2, second odrer effects can be ignored if the 
% slenderness ratio is lower than the limiting slenderness.
%
% Input:    Information about the member (structure array)
% Output:   A logical value (0 or 1)
%--------------------------------------------------------------------------

% {Analysis of second order effects with axial load, General - 5.8.2}
% (Needs to be read)

% {5.8.2 (6)}
% If second ordre effects are less than 10% of the first order effects,
% then they can be ignored.

% {Slenderness criterion for isolated members - 5.8.3.1 (1) - p.65}

    % Intermediate variables
    A = 1 / (1 + 0.2*M.varphi_ef);
    B = sqrt(1 + 2*M.omega);
    C = 1.7 - M.r_m;

    % The limiting slenderness 
    lambda_lim = 20*A*B*C*sqrt( ( M.A_c*M.f_cd )/M.N_Ed );

% Comparing to the column's slenderness ratios
if M.lambda < lambda_lim

    % If the limiting slenderness is greater than the actual slenderness,
    % a value is assigned accordingly.
    soe = 0;

else

    % The procedure is the same, if the opposite is true 
    soe = 1;

end


%--------------------------------------------------------------------------
%% Outputting requested information


if info_lev == 1
    fprintf('(lambda < lambda_lim) = (%.2f < %.2f) [-] | ' ...
        , M.lambda, lambda_lim)

elseif info_lev == 2
    
    % Information depends on the result
    if soe == 1
        fprintf(['The second order effects can NOT neglected, ' ...
            'becasuse: (lambda > lambda_lim) = (%.2f > %.2f) [-] \n'] ...
            , M.lambda, lambda_lim)
    else
        fprintf(['The second order effects CAN neglected, ' ...
            'becasuse: (lambda < lambda_lim) = (%.2f < %.2f) [-] \n'] ...
            , M.lambda, lambda_lim)
    end

end