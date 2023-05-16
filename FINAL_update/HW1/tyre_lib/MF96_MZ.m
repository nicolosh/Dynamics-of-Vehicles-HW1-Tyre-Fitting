% Self Aligning Moment
function [MZ] = MF96_MZ(kappa, alpha, phi, Fz, tyre_data)

 % precode

  [FY0] = MF96_FY0(kappa, alpha, phi, Fz, tyre_data);
  [MZr] = MF96_MZr(kappa, alpha, phi, Fz, tyre_data);
  [t] = MF96_t(kappa, alpha, phi, Fz, tyre_data);

 % main code

  MZ = -FY0 * t + MZr;
  
 end
