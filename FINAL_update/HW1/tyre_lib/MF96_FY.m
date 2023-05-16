% Combined lateral force
function [FY] = MF96_FY(kappa, alpha, phi, Fz, tyre_data)

 % precode

  [Gxa, Gyk, SVyk] = MF96_FXFYCOMB_coeffs(kappa, alpha, phi, Fz, tyre_data);
  [Fy0] = MF96_FY0(kappa, alpha, phi, Fz, tyre_data);

 % main code

  FY = Gyk * Fy0 + SVyk;
  
 end
