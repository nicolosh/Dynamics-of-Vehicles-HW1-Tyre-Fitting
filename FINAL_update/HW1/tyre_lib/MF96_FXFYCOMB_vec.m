% Pure lateral force FY0
% this function remap the scalar function to its vectorial form
function [Gxa_vec,Gyk_vec,SVyk_vec] = MF96_FXFYCOMB_vec(kappa_vec, alpha_vec, phi_vec, Fz_vec, tyre_data)

  Gxa_vec = zeros(size(alpha_vec));
  Gyk_vec = zeros(size(alpha_vec));
  SVyk_vec = zeros(size(alpha_vec));
  for i = 1:length(alpha_vec)
   
     [Gxa_vec(i),Gyk_vec(i),SVyk_vec(i)] = MF96_FXFYCOMB_coeffs(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);
   
  end
  
 end
