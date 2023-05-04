% Combined longitudinal force
function [FX] = MF96_FX(kappa, alpha, phi, Fz, tyre_data)

 % precode

  [Gxa, Gyk, SVyk] = MF96_FXFYCOMB_coeffs(kappa, alpha, phi, Fz, tyre_data);
  [Fx0] = MF96_FX0(kappa, alpha, phi, Fz, tyre_data);

 % main code

  FX = Gxa * Fx0;
  
 end
