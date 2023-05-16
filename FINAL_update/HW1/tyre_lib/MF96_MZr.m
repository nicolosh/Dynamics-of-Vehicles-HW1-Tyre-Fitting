% Residual Torque
function [MZr] = MF96_MZr(kappa, alpha, phi, Fz, tyre_data)

 % precode

  [alpha__r, Br, Dr, Bt, Ct, Dt, Et, alpha__t] = MF96_MZ0_coeffs(kappa, alpha, phi, Fz, tyre_data);

 % main code

  MZr = Dr * ((Br ^ 2 * alpha__r ^ 2 + 1) ^ (-0.1e1 / 0.2e1)) * cos(alpha);
  
 end
