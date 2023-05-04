function res = resid_pure_Mz_varGamma(P, MZ, ALPHA, GAMMA, FZ, FY_vec, tyre_data)

    % ----------------------------------------------------------------------
    %% Compute the residuals - least squares approach - to fit the Mz curve 
    %  with Fz=Fz_nom, IA=0. Pacejka 1996 Magic Formula
    % ----------------------------------------------------------------------

    % Define MF coefficients

    %Fz0 = 200*4.44822; % Nominal load 200 lbf
    tmp_tyre_data = tyre_data;
	
    % assign computed parameter
   
    tmp_tyre_data.qBz4 = P(1); 
    tmp_tyre_data.qBz5 = P(2); 
    tmp_tyre_data.qEz5 = P(3);  
        
    % Longitudinal Force (Pure Longitudinal Slip) Equations
    res = 0;
    for i=1:length(ALPHA)
       mz0  = MF96_MZ0(ALPHA(i), GAMMA(i), FZ, FY_vec(i), tmp_tyre_data);
       res = res+(mz0-MZ(i))^2;
    end
    
    % Compute the residuals
    res = res/sum(MZ.^2);

end

