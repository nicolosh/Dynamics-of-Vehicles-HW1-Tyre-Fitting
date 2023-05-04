% Coefficients for Magic Formula lateral/longitudinal combined forces
function [Gxa, Gyk, SVyk] = MF96_FXFYCOMB_coeffs(kappa, alpha, phi, Fz, tyre_data)

 % precode

  rBx1            = tyre_data.rBx1;
  rBx2            = tyre_data.rBx2;
  rBy1            = tyre_data.rBy1;
  rBy2            = tyre_data.rBy2;
  rBy3            = tyre_data.rBy3;
  rCx1            = tyre_data.rCx1;
  rCy1            = tyre_data.rCy1;
  rHx1            = tyre_data.rHx1;
  rHy1            = tyre_data.rHy1;
  rVy1            = tyre_data.rVy1;
  rVy2            = tyre_data.rVy2;
  rVy3            = tyre_data.rVy3;
  rVy4            = tyre_data.rVy4;
  rVy5            = tyre_data.rVy5;
  rVy6            = tyre_data.rVy6;
  LFZ0            = tyre_data.LFZ0;
  LGAMMAY         = tyre_data.LGAMMAY;
  LMUY            = tyre_data.LMUY;
  LVYK            = tyre_data.LVYK;
  LXA             = tyre_data.LXA;
  LYK             = tyre_data.LYK;
  FZ0             = tyre_data.FZ0;
  Fz01            = tyre_data.Fz01;
  pCy1            = tyre_data.pCy1;
  pDy1            = tyre_data.pDy1;
  pDy2            = tyre_data.pDy2;
  pDy3            = tyre_data.pDy3;
  pEy1            = tyre_data.pEy1;
  pEy2            = tyre_data.pEy2;
  pEy3            = tyre_data.pEy3;
  pEy4            = tyre_data.pEy4;
  pHy1            = tyre_data.pHy1;
  pHy2            = tyre_data.pHy2;
  pHy3            = tyre_data.pHy3;
  pKy1            = tyre_data.pKy1;
  pKy2            = tyre_data.pKy2;
  pKy3            = tyre_data.pKy3;
  pVy1            = tyre_data.pVy1;
  pVy2            = tyre_data.pVy2;
  pVy3            = tyre_data.pVy3;
  pVy4            = tyre_data.pVy4;
  

 % main code

  FZ01 = (LFZ0 * FZ0);
  dfz = Fz / FZ01 - 1;
  SHxa = rHx1;
  Bxa = rBx1 * (kappa ^ 2 * rBx2 ^ 2 + 1) ^ (-0.1e1 / 0.2e1) * LXA;
  Cxa = rCx1;
  Dxa = 0.1e1 / cos(Cxa * atan((Bxa * SHxa)));
  Gxa = Dxa * cos(Cxa * atan((Bxa * (alpha + SHxa))));
  gamma__s = (phi * LGAMMAY);
  mu__y = (dfz * pDy2 + pDy1) * (-pDy3 * gamma__s ^ 2 + 1) * LMUY;
  SHyk = rHy1;
  DVyk = mu__y * Fz * (dfz * rVy2 + phi * rVy3 + rVy1) * (alpha ^ 2 * rVy4 ^ 2 + 1) ^ (-0.1e1 / 0.2e1);
  SVyk = DVyk * sin(rVy5 * atan((rVy6 * kappa))) * LVYK;
  Byk = rBy1 * (1 + rBy2 ^ 2 * (alpha - rBy3) ^ 2) ^ (-0.1e1 / 0.2e1) * LYK;
  Cyk = rCy1;
  Dyk = 0.1e1 / cos(Cyk * atan((Byk * SHyk)));
  Gyk = Dyk * cos(Cyk * atan((Byk * (kappa + SHyk))));
  
 end
