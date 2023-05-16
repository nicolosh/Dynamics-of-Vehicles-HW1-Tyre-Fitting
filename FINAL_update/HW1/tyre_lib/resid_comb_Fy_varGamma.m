function res = resid_comb_Fy_dfz(P,FY,KAPPA,ALPHA,GAMMA,FZ,tyre_data)

    % ----------------------------------------------------------------------
    %% Compute the residuals - least squares approach - to fit the Fy curve 
    %  with Fz=Fz_nom, IA=0. Pacejka 1996 Magic Formula
    % ----------------------------------------------------------------------

    % Define MF coefficients

    %Fz0 = 200*4.44822; % Nominal load 200 lbf
    
    tmp_tyre_data = tyre_data;
    
       
    tmp_tyre_data.rVy3 = P(1); 
    
    
   %dfz = (Z - Fz0)./Fz0 ;
    
    % Longitudinal Force (Pure Longitudinal Slip) Equations
    res = 0;
    for i=1:length(ALPHA)
       fy  = MF96_FY(KAPPA(i), ALPHA(i), GAMMA(i), FZ(i), tmp_tyre_data);
       res = res+(fy-FY(i))^2;
    end
    
    % Compute the residuals
    res = res/sum(FY.^2);

end
