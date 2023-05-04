% Pure lateral force FY0
% this function remap the scalar function to its vectorial form
function [mz0_vec] = MF96_MZ0_vec(kappa_vec, alpha_vec, phi_vec, Fz_vec, Fy_vec, tyre_data)

  
  mz0_vec = zeros(size(alpha_vec));
  for i = 1:length(alpha_vec)
   % precode
   [alpha__r, Br, Dr, Bt, Ct, Dt, Et, alpha__t] = MF96_MZ0_coeffs(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);
   
   %[FY] = MF96_FY(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);
   
   [t] = MF96_t(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);
   
   [MZr] = MF96_MZr(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);

   % main code
    mz0_vec(i) = -FY_vec(i) * t + MZr;
  end
  
 end
