% Self Aligning Moment
function [MZ0] = MF96_MZ0(kappa, alpha, phi, Fz, tyre_data)

 % precode

  [alpha__r, Br, Dr, Bt, Ct, Dt, Et, alpha__t] = MF96_MZ0_coeffs(kappa, alpha, phi, Fz, tyre_data);
  [FY] = MF96_FY(kappa, alpha, phi, Fz, tyre_data);
  [t] = MF96_t(kappa, alpha, phi, Fz, tyre_data);
  [MZr] = MF96_MZr(kappa, alpha, phi, Fz, tyre_data);

 % main code

  MZ0 = -Fy0 * t + MZr;
  
 end
