function varphi_ef = effective_creep(info_lev, M)

%--------------------------------------------------------------------------
% Determines the effective creep coefficient for outside conditions. 
% Firstly a check is conducted, in accordance with EC2, to see if the 
% effects of creep can be ignored.
%
% Assumptions:
%  - It cement class N is assumed, because it results in the highest final
%    creep coefficient of all cements produced by AalborgPortland.
%  - Temperature changes from casting to time of loading is negligible.
%  - Every surface is subject to drying, and not obstruct.
% 
% {Creep - 5.8.4 (3)}
% The ratio is calculated in the section with the maximum moment, if 
% M_0Eqp/M_0Ed varies through the column.
%
% {Annex B}
%
% Input:    Information about the member (structure array)
% Output:   Effective creep coefficient
%--------------------------------------------------------------------------

% Relative humidity (RH) for outside conditions
RH = 80;                                            %[%]    {3.1.4 (5)}

% Preliminary calculation 

    % Factor for cement strength class
    alpha = 0;                                      %[-]    {B.1 (B.9)}

    % The notional size of the member, theoretical size
    h_0 = (2*M.A_c)/(2*M.h + 2*M.w) *10^(3);        %[mm]   {3.1.4 (5)}

    % Mean value of compressive concrete strength
    f_cm = M.f_ck*10^(-6) + 8;                      %[MPa]  {3.1.3 (T)}


% Reduction and correction factors

    % Correcting the concrete's age at loading
    t_0 = max(M.t0*(9/(2+M.t0^(1.2)))^alpha ,0.5);  %[days] {B.1 (B.9)}

    % Factor consider effects of RH on creep coe.
    varphi_RH = 1 + (1-RH/100)/(0.1*h_0^(1/3));     %[-]    {B.1 (B.3a)}

    % Correction factor for concrete strength on notional creep coe.
    beta_fcm = 16.8/sqrt(f_cm);                     %[-]    {B.1 (B.4)}

    % Correction factor for concrete age on notional creep coe.
    beta_t0 = 1/(0.1 + t_0^(0.20));                 %[-]    {B.1 (B.5)}

    % Coefficient describing the development of creep with time, also
    % considering the concrete age at loading, if t = infty.
    beta_c = 1;

    
    % Changing factors if the mean strength is greater than 35
    if f_cm > 35
        
        % Reduction factors                         %[-]    {B.1 (B.8c)}
        alpha_1 = (35/f_cm)^(0.7);
        alpha_2 = (35/f_cm)^(0.2);

        % Effects of RH on creep coe.               %[-]    {B.1 (B.3b)}
        varphi_RH = (1 + (1-RH/100)/(0.1*h_0^(1/3)) * alpha_1) * alpha_2;

    end


% Final calcualtions

    % Notional creep coefficient
    varphi_0 = varphi_RH * beta_fcm * beta_t0;      %[-]        {B.1 (B.2)}

    % Creep coefficent
    varphi = varphi_0 * beta_c;                     %[-]        {B.1 (B.1)}


% Checks the three conditions, if all are satisfied creep can be ignored.
if (varphi <= 2) && (M.lambda <= 75) && (M.M_0Ed/M.N_Ed >= M.h)

    % Sets the effective creep coefficient to zero.
    varphi_ef = 0;                                          %{5.8.4 (4)}

else

    % If one or more condition is not satisfied, then the effective
    % creep coefficient is calculated in accordance with:
    % {Creep - 5.8.4 (2)}
    varphi_ef = varphi * M.M_0Eqp/M.M_0Ed;

end

%--------------------------------------------------------------------------
%% Outputting requested information

if info_lev == 1
    fprintf(' varphi_ef = %.2e [-] | varphi_0 = %.2e [-] | ', ...
        varphi_ef, varphi_0)

elseif info_lev == 2
    
    % Information depends on the result
    if varphi_ef == 0
        fprintf(['The effective creep coefficient CAN be neglected ' ...
            'and is therefore set to: %.2e [-] \n'], varphi_ef);
        fprintf(['Becasue: \n' ...
            ' - varphi = %.2f <= 2 [-] \n' ...
            ' - lambda = %.2f <= 75 [-] \n' ...
            ' - M_0Ed/N_Ed = %.2f/%.2f = %.2f >= h = %.2f [m]\n'] ...
            , varphi, M.lambda, M.M_0Ed, M.N_Ed ...
            , M.M_0Ed/M.N_Ed, M.h)

    else
        fprintf(['The effective creep coefficient can NOT be ' ...
            'neglected and is therefore set to: %.2e [-] \n'], varphi_ef);
        fprintf('Bacause: \n');

        if (varphi > 2)
            fprintf(' - varphi = %.2f > 2 [-] \n', varphi)
        end

        if (M.lambda > 75)
            fprintf(' - lambda = %.2f > 2 [-] \n', M.lambda)
        end

        if (M.M_0Ed/M.N_Ed < M.h)
            fprintf([' - M_0Ed/N_Ed = %.f/%.f = %.4f < ' ...
                'h = %.4f [m] \n'], M.M_0Ed, M.N_Ed ...
            , M.M_0Ed/M.N_Ed, M.h)
        end

    end

end