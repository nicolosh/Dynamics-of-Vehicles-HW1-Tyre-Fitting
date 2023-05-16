% Residual Torque
function [t] = MF96_t(kappa, alpha, phi, Fz, tyre_data)

 % precode

  [alpha__r, Br, Dr, Bt, Ct, Dt, Et, alpha__t] = MF96_MZ0_coeffs(kappa, alpha, phi, Fz, tyre_data);

 % main code

  t = Dt * cos(Ct * atan(-Bt * alpha__t + Et * (Bt * alpha__t - atan(Bt * alpha__t)))) * cos(alpha);
  
 end
