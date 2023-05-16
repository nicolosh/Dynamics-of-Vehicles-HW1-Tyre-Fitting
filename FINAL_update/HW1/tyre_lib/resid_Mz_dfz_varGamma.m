function res = resid_Mz_dfz_varGamma(P,MZ,ALPHA,GAMMA,FZ,tyre_data)

    % ----------------------------------------------------------------------
    %% Compute the residuals - least squares approach - to fit the Mz curve 
    %  with Fz=var, IA=0. Pacejka 1996 Magic Formula
    % ----------------------------------------------------------------------

    % Define MF coefficients

    %Fz0 = 200*4.44822; % Nominal load 200 lbf
    
    tmp_tyre_data = tyre_data;
   
    tmp_tyre_data.qHz4 = P(1);
    tmp_tyre_data.qDz9 = P(2);
    

    
    % Longitudinal Force (Pure Longitudinal Slip) Equations
    res = 0;
    for i=1:length(ALPHA)
       Mz0  = MF96_MZ(0, ALPHA(i), GAMMA(i), FZ(i), tmp_tyre_data);
       res = res+(Mz0-MZ(i))^2;
       %res = res+(fx0/FX(i)-1)^2;
    end
    
    % Compute the residuals
    res = res/sum(MZ.^2);
    

end