% Self Aligning moment
% this function remap the scalar function to its vectorial form
function [mz_vec] = MF96_MZ_vec(kappa_vec, alpha_vec, phi_vec, Fz_vec, tyre_data)

  
  mz_vec = zeros(size(alpha_vec));
  for i = 1:length(alpha_vec)
         mz_vec(i) = MF96_MZ(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i),tyre_data);
  end
  
 end
