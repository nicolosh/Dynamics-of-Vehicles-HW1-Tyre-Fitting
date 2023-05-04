%% Initialisation
clc
clearvars 
close all   

% Set LaTeX as default interpreter for axis labels, ticks and legends
set(0,'defaulttextinterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

set(0,'DefaultFigureWindowStyle','docked');
set(0,'defaultAxesFontSize',  16)
set(0,'DefaultLegendFontSize',16)

addpath('tyre_lib/')


to_rad = pi/180;
to_deg = 180/pi;

%% LATERAL
%% Select tyre dataset
%dataset path
data_set_path = 'dataset/';
% dataset selection and loading

data_set = 'Hoosier_B1464run23'; % pure lateral forces
%data_set = 'Hoosier_B1464run30';  % braking/traction (pure log. force) + combined

% tyre geometric data:
% Hoosier	18.0x6.0-10
% 18 diameter in inches
% 6.0 section width in inches
% tread width in inches
diameter = 18*2.56; %
Fz0 = 220;   % [N] nominal load is given
R0  = diameter/2/100; % [m] get from nominal load R0 (m) *** TO BE CHANGED ***


fprintf('Loading dataset ...')
switch data_set
    case 'Hoosier_B1464run23'
  load ([data_set_path, 'Hoosier_B1464run23.mat']); % pure lateral
  cut_start = 31300;
  cut_end   = 54500;
    case 'Hoosier_B1464run30'
  load ([data_set_path, 'Hoosier_B1464run30.mat']); % pure longitudinal
  cut_start = 19028;
  cut_end   = 37643;
  otherwise 
  error('Not found dataset: `%s`\n', data_set) ;
  
end

% select dataset portion
smpl_range = cut_start:cut_end;

fprintf('completed!\n')
%% Plot raw data

figure
tiledlayout(6,1)

ax_list(1) = nexttile; y_range = [min(min(-FZ),0) round(max(-FZ)*1.1)];
plot(-FZ)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')

ax_list(2) = nexttile; y_range = [min(min(IA),0) round(max(IA)*1.1)];
plot(IA)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(3) = nexttile; y_range = [min(min(SA),0) round(max(SA)*1.1)];
plot(SA)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Side slip')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(4) = nexttile; y_range = [min(min(SL),0) round(max(SL)*1.1)];
plot(SL)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Longitudinal slip')
xlabel('Samples [-]')
ylabel('[-]')

ax_list(5) = nexttile; y_range = [min(min(P),0) round(max(P)*1.1)];
plot(P)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Tyre pressure')
xlabel('Samples [-]')
ylabel('[psi]')

ax_list(6) = nexttile;  y_range = [min(min(TSTC),0) round(max(TSTC)*1.1)];
plot(TSTC,'DisplayName','Center')
hold on
plot(TSTI,'DisplayName','Internal')
plot(TSTO,'DisplayName','Outboard')
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Tyre temperatures')
xlabel('Samples [-]')
ylabel('[degC]')

linkaxes(ax_list,'x')


%plot(SA,FY)

%% Select some specific data
% Cut crappy data and select only 12 psi data
clear idx;
vec_samples = 1:1:length(smpl_range);
tyre_data = table(); % create empty table
% store raw data in table
tyre_data.SL =  SL(smpl_range);
tyre_data.SA =  SA(smpl_range)*to_rad;
tyre_data.FZ = -FZ(smpl_range);  % 0.453592  lb/kg
tyre_data.FX =  FX(smpl_range);
tyre_data.FY =  FY(smpl_range);
tyre_data.MZ =  MZ(smpl_range);
tyre_data.IA =  IA(smpl_range)*to_rad;

% Extract points at constant inclination angle
GAMMA_tol = 0.05*to_rad;
idx.GAMMA_0 = 0.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 0.0*to_rad+GAMMA_tol;
idx.GAMMA_1 = 1.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 1.0*to_rad+GAMMA_tol;
idx.GAMMA_2 = 2.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 2.0*to_rad+GAMMA_tol;
idx.GAMMA_3 = 3.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 3.0*to_rad+GAMMA_tol;
idx.GAMMA_4 = 4.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 4.0*to_rad+GAMMA_tol;
idx.GAMMA_5 = 5.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 5.0*to_rad+GAMMA_tol;
GAMMA_0  = tyre_data( idx.GAMMA_0, : );
GAMMA_1  = tyre_data( idx.GAMMA_1, : );
GAMMA_2  = tyre_data( idx.GAMMA_2, : );
GAMMA_3  = tyre_data( idx.GAMMA_3, : );
GAMMA_4  = tyre_data( idx.GAMMA_4, : );
GAMMA_5  = tyre_data( idx.GAMMA_5, : );

% Extract points at constant vertical load
% Test data done at: 
%  - 50lbf  ( 50*0.453592*9.81 =  223N )
%  - 150lbf (150*0.453592*9.81 =  667N )
%  - 200lbf (200*0.453592*9.81 =  890N )
%  - 250lbf (250*0.453592*9.81 = 1120N )

FZ_tol = 100;
idx.FZ_220  = 220-FZ_tol < tyre_data.FZ & tyre_data.FZ < 220+FZ_tol;
idx.FZ_440  = 440-FZ_tol < tyre_data.FZ & tyre_data.FZ < 440+FZ_tol;
idx.FZ_700  = 700-FZ_tol < tyre_data.FZ & tyre_data.FZ < 700+FZ_tol;
idx.FZ_900  = 900-FZ_tol < tyre_data.FZ & tyre_data.FZ < 900+FZ_tol;
idx.FZ_1120 = 1120-FZ_tol < tyre_data.FZ & tyre_data.FZ < 1120+FZ_tol;
FZ_220  = tyre_data( idx.FZ_220, : );
FZ_440  = tyre_data( idx.FZ_440, : );
FZ_700  = tyre_data( idx.FZ_700, : );
FZ_900  = tyre_data( idx.FZ_900, : );
FZ_1120 = tyre_data( idx.FZ_1120, : );


figure()
tiledlayout(3,1)

ax_list(1) = nexttile;
plot(tyre_data.IA*to_deg)
hold on
plot(vec_samples(idx.GAMMA_0),GAMMA_0.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_1),GAMMA_1.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_2),GAMMA_2.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_3),GAMMA_3.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_4),GAMMA_4.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_5),GAMMA_5.IA*to_deg,'.');
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(2) = nexttile;
plot(tyre_data.FZ)
hold on
plot(vec_samples(idx.FZ_220),FZ_220.FZ,'.');
plot(vec_samples(idx.FZ_440),FZ_440.FZ,'.');
plot(vec_samples(idx.FZ_700),FZ_700.FZ,'.');
plot(vec_samples(idx.FZ_900),FZ_900.FZ,'.');
plot(vec_samples(idx.FZ_1120),FZ_1120.FZ,'.');
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')


ax_list(3) = nexttile;
plot(tyre_data.SA)
hold on
title('Slide slip')
xlabel('Samples [-]')
ylabel('[rad]')

%% -------------------------------------------------------------------------
% FITTING WITH GUESS VALUES and nominal vertical load
%--------------------------------------------------------------------------
%% Intersect tables to obtain specific sub-datasets

[TData1, ~] = intersect_table_data( GAMMA_0, FZ_220 );

%% FITTING 
% initialise tyre data
tyre_coeffs = initialise_tyre_data(R0, Fz0);

%% Fitting with Fz=Fz_nom= 220N and camber=0 VX= 10
% ------------------
% lateral slip
clear P;
% Fit the coeffs {pCy1, pDy1, pEy1, pEy4, pKy1, pHy1, pVy1}
FZ0 = mean(TData1.FZ);

zeros_vec = zeros(size(TData1.SA));
ones_vec  = ones(size(TData1.SA));

% k = 0, alpha, phi = 0 (camber)
FY0_guess = MF96_FY0_vec(zeros_vec, TData1.SA, zeros_vec, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);

% check guess 
figure('Name','FY0 guess')
plot(TData1.SA,TData1.FY,'.')
hold on
plot(TData1.SA,FY0_guess,'-')

% Plot raw data and initial guess
% figure()
% plot(TDataSub.KAPPA,TDataSub.FX,'o')
% hold on
% plot(TDataSub.KAPPA,FX0_guess,'x')

%% ------------------------------------------------------------------------
% FITTING WITH NOMINAL LOAD Fz=Fz_nom = 220N and camber = 0  alpha = 0 VX = 10
%--------------------------------------------------------------------------
% Guess values for parameters to be optimised
%    [pCy1 pDy1 pEy1 pHy1  pKy1  pKy2  pVy1] 
P0 = [0.1,0.1,0.1,0.1,0.1,0.1,0.1]; 


% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
%    [pCy1 pDy1 pEy1 pHy1  pKy1  pKy2  pVy1]
lb = [1,-1000,-1000,-1000,-1000,-1000,-1000];
ub = [2,1000,1000,1000,1000,1000,1000];

ALPHA_vec = TData1.SA; % side slip angle vector
FY_vec = TData1.FY; % lateral force

% check guess

% FY0_fz_nom_vec = MF96_FY0_vec(SA_vec,zeros(size(SA_vec)) , zeros(size(SA_vec)), ...
%                               FZ0.*ones(size(SA_vec)),tyre_coeffs);

% residuals

[P_fz_nom,fval,exitflag,iterations] = fmincon(@(P)resid_pure_Fy(P,FY_vec, ALPHA_vec,0,FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);


% Update tyre data with new optimal values                             
tyre_coeffs.pCy1 = P_fz_nom(1) ; % 1
tyre_coeffs.pDy1 = P_fz_nom(2) ;  
tyre_coeffs.pEy1 = P_fz_nom(3) ;
tyre_coeffs.pHy1 = P_fz_nom(4) ;
tyre_coeffs.pKy1 = P_fz_nom(5) ; 
tyre_coeffs.pKy2 = P_fz_nom(6) ;
tyre_coeffs.pVy1 = P_fz_nom(7) ;

% lateral side slip to be used in the plot for checking the interpolation outside the input data range
SA_vec = -0.25:0.001:0.25; 

% zero camber, zero k and alpha
FY0_fz_nom_vec = MF96_FY0_vec(zeros(size(SA_vec)), SA_vec, zeros(size(SA_vec)), ...
                              FZ0.*ones(size(SA_vec)),tyre_coeffs);


figure('Name','Fy0(Fz0)')
plot(TData1.SA,TData1.FY,'o')
hold on
plot(SA_vec,FY0_fz_nom_vec,'-','LineWidth',2)
xlabel('$\alpha$ [rad]')
ylabel('$F_{y0}$ [N]')

%% ------------------------------------------------------------------------
% FIT COEFFICIENTS WITH VARIABLE LOAD
%--------------------------------------------------------------------------
%% Fit coefficient with variable load
% extract data with variable load
% [TDataDFz1, ~] = intersect_table_data( GAMMA_0 );

% Lateral slip

% Fit the coeffs {pCx1, pDx1, pEx1, pEx4, pKx1, pHx1, pVx1}
%FZ0 = mean(TData0.FZ);

zeros_vec = zeros(size(GAMMA_0.SA));
ones_vec  = ones(size(GAMMA_0.SA));

Fy0_guess = MF96_FY0_vec(zeros_vec, GAMMA_0.SA, zeros_vec, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);

% check guess 
figure()
plot(GAMMA_0.SA,GAMMA_0.FY,'.')
hold on
plot(GAMMA_0.SA,Fy0_guess,'-')

% Plot raw data and initial guess
% figure()
% plot(TDataSub.KAPPA,TDataSub.FX,'o')
% hold on
% plot(TDataSub.KAPPA,FX0_guess,'x')

% Guess values for parameters to be optimised
%    [pDx2 pEx2 pEx3 pHx2  pKx2  pKx3  pVx2] 
P0 = [0.1, 0.1, 0.1, 0.1]; 


% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
% 1< pCx1 < 2 
% 0< pEx1 < 1 
%    [pDy2 pEy2 pHy2 pVy2] 
lb = [-1000,0,-1000,-1000];
ub = [1000,1,1000,1000];


ALPHA_vec = GAMMA_0.SA;
FY_vec    = GAMMA_0.FY;
FZ_vec    = GAMMA_0.FZ;

% check guess
% SL_vec = -0.3:0.001:0.3;
% FX0_dfz_vec = MF96_FX0_vec(SL_vec,zeros(size(SL_vec)) , zeros(size(SL_vec)), ...
%                            TDataDFz.FZ,tyre_coeffs);
% 
% figure
% plot(KAPPA_vec,FX_vec,'.')
% hold on
% plot(SL_vec,FX0_dfz_vec,'.')


% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 

[P_dfz,fval,exitflag] = fmincon(@(P)resid_pure_Fy_varFz(P,FY_vec, ALPHA_vec,0,FZ_vec, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

disp(exitflag)
% Change tyre data with new optimal values                             
tyre_coeffs.pDy2 = P_dfz(1) ; % 1
tyre_coeffs.pEy2 = P_dfz(2) ;  
tyre_coeffs.pHy2 = P_dfz(3) ;
tyre_coeffs.pVy2 = P_dfz(4) ;


res_FY0_dfz_vec = resid_pure_Fy_varFz(P_dfz,FY_vec,SA_vec,0 , FZ_vec,tyre_coeffs);

tmp_zeros = zeros(size(SA_vec));
tmp_ones = ones(size(SA_vec));


FY0_fz_var_vec1 = MF96_FY0_vec(tmp_zeros, SA_vec ,tmp_zeros, mean(FZ_220.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec2 = MF96_FY0_vec(tmp_zeros, SA_vec ,tmp_zeros, mean(FZ_440.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec3 = MF96_FY0_vec(tmp_zeros, SA_vec ,tmp_zeros, mean(FZ_700.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec4 = MF96_FY0_vec(tmp_zeros, SA_vec ,tmp_zeros, mean(FZ_900.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec5 = MF96_FY0_vec(tmp_zeros, SA_vec ,tmp_zeros, mean(FZ_1120.FZ)*tmp_ones,tyre_coeffs);


figure('Name','Fy0(Fz0)')
plot(GAMMA_0.SA,GAMMA_0.FY,'o')
hold on
%plot(TDataSub.KAPPA,FX0_fz_nom_vec,'-')
%plot(SL_vec,FX0_dfz_vec,'-','LineWidth',2)
plot(SA_vec,FY0_fz_var_vec1,'-','LineWidth',2)
plot(SA_vec,FY0_fz_var_vec2,'-','LineWidth',2)
plot(SA_vec,FY0_fz_var_vec3,'-','LineWidth',2)
plot(SA_vec,FY0_fz_var_vec4,'-','LineWidth',2)
plot(SA_vec,FY0_fz_var_vec5,'-','LineWidth',2)

xlabel('$\alpha$ [rad]')
ylabel('$F_{y0}$ [N]')

% 
% [kappa__x, Bx, Cx, Dx, Ex, SVx] =MF96_FX0_coeffs(0, 0, 0, mean(FZ_220.FZ), tyre_coeffs);
% Calfa_vec1_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);
% [kappa__x, Bx, Cx, Dx, Ex, SVx] =MF96_FX0_coeffs(0, 0, 0, mean(FZ_700.FZ), tyre_coeffs);
% Calfa_vec2_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);
% [kappa__x, Bx, Cx, Dx, Ex, SVx] =MF96_FX0_coeffs(0, 0, 0, mean(FZ_900.FZ), tyre_coeffs);
% Calfa_vec3_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);
% [kappa__x, Bx, Cx, Dx, Ex, SVx] =MF96_FX0_coeffs(0, 0, 0, mean(FZ_1120.FZ), tyre_coeffs);
% Calfa_vec4_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);
% 
% Calfa_vec1 = MF96_CorneringStiffness(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_220.FZ)*tmp_ones,tyre_coeffs);
% Calfa_vec2 = MF96_CorneringStiffness(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_700.FZ)*tmp_ones,tyre_coeffs);
% Calfa_vec3 = MF96_CorneringStiffness(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_900.FZ)*tmp_ones,tyre_coeffs);
% Calfa_vec4 = MF96_CorneringStiffness(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_1120.FZ)*tmp_ones,tyre_coeffs);
% 
% figure('Name','C_alpha')
% subplot(2,1,1)
% hold on
% %plot(TDataSub.KAPPA,FX0_fz_nom_vec,'-')
% plot(mean(FZ_220.FZ),Calfa_vec1_0,'+','LineWidth',2)
% plot(mean(FZ_700.FZ),Calfa_vec3_0,'+','LineWidth',2)
% plot(mean(FZ_900.FZ),Calfa_vec4_0,'+','LineWidth',2)
% plot(mean(FZ_1120.FZ),Calfa_vec2_0,'+','LineWidth',2)
% legend({'$Fz_{220}$','$Fz_{700}$','$Fz_{900}$','$Fz_{1120}$'})
% 
% subplot(2,1,2)
% hold on
% %plot(TDataSub.KAPPA,FX0_fz_nom_vec,'-')
% plot(SL_vec,Calfa_vec1,'-','LineWidth',2)
% plot(SL_vec,Calfa_vec2,'-','LineWidth',2)
% plot(SL_vec,Calfa_vec3,'-','LineWidth',2)
% plot(SL_vec,Calfa_vec4,'-','LineWidth',2)
% legend({'$Fz_{220}$','$Fz_{700}$','$Fz_{900}$','$Fz_{1120}$'})
% 

%% ------------------------------------------------------------------------
% FIT COEFFICIENTS WITH VARIABLE CAMBER
%--------------------------------------------------------------------------
%% Fit coefficient with variable camber

% extract data with variable load
% [TDataGamma1, ~] = intersect_table_data( FZ_220 );

% Fit the coeffs { pDx3}

% Guess values for parameters to be optimised
P0 = [0.1,0.1,0.1,0.1,0.1,0.1,0.1]; 

% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
% 1< pCx1 < 2 
% 0< pEx1 < 1 
%lb = [0, 0,  0, 0,  0,  0,  0];
%ub = [2, 1e6,1, 1,1e1,1e2,1e2];
lb = [];
ub = [];


zeros_vec = zeros(size(FZ_220.SA));
ones_vec  = ones(size(FZ_220.SA));

ALPHA_vec = FZ_220.SA;
GAMMA_vec = FZ_220.IA; 
FY_vec    = FZ_220.FY;
FZ_vec    = FZ_220.FZ;

figure()
plot(ALPHA_vec,FY_vec);


% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
[P_varGamma1,fval,exitflag] = fmincon(@(P)resid_pure_Fy_varGamma(P,FY_vec, ALPHA_vec,GAMMA_vec,tyre_coeffs.FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Change tyre data with new optimal values                             
tyre_coeffs.pDy3 = P_varGamma1(1) ; % 1
tyre_coeffs.pEy3 = P_varGamma1(2) ;
tyre_coeffs.pEy4 = P_varGamma1(3) ;
tyre_coeffs.pHy3 = P_varGamma1(4) ;
tyre_coeffs.pKy3 = P_varGamma1(5) ;
tyre_coeffs.pVy3 = P_varGamma1(6) ;
tyre_coeffs.pVy4 = P_varGamma1(7) ;



% FY0_varGamma_vec = MF96_FY0_vec(zeros_vec, SA_vec, GAMMA_vec, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);

FY0_varGamma_vec0 = MF96_FY0_vec(zeros_vec, SA_vec, GAMMA_0.IA, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);
FY0_varGamma_vec1 = MF96_FY0_vec(zeros_vec, SA_vec, GAMMA_1.IA, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);
FY0_varGamma_vec2 = MF96_FY0_vec(zeros_vec, SA_vec, GAMMA_2.IA, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);
FY0_varGamma_vec3 = MF96_FY0_vec(zeros_vec, SA_vec, GAMMA_3.IA, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);
FY0_varGamma_vec4 = MF96_FY0_vec(zeros_vec, SA_vec, GAMMA_4.IA, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);

figure('Name','Fy0 vs Gamma')
plot(ALPHA_vec,FZ_220.FY,'o')
hold on
plot(SA_vec,FY0_varGamma_vec0,'-')
plot(SA_vec,FY0_varGamma_vec1,'-')
plot(SA_vec,FY0_varGamma_vec2,'-')
plot(SA_vec,FY0_varGamma_vec3,'-')
plot(SA_vec,FY0_varGamma_vec4,'-')
xlabel('$\alpha$ [rad]')
ylabel('$F_{y0}$ [N]')
% Calculate the residuals with the optimal solution found above
%res_Fx0_varGamma  = resid_pure_Fx_varGamma(P_varGamma,FX_vec, KAPPA_vec,GAMMA_vec,tyre_coeffs.FZ0, tyre_coeffs);

% R-squared is 
% 1-SSE/SST
% SSE/SST = res_Fx0_nom

% % SSE is the sum of squared error,  SST is the sum of squared total
% fprintf('R-squared = %6.3f\n',1-res_Fx0_varGamma);
% 
% 
% [kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(0, 0, GAMMA_vec(3), tyre_coeffs.FZ0, tyre_coeffs);
% % 
% fprintf('Bx      = %6.3f\n',Bx);
% fprintf('Cx      = %6.3f\n',Cx);
% fprintf('mux      = %6.3f\n',Dx/tyre_coeffs.FZ0);
% fprintf('Ex      = %6.3f\n',Ex);
% fprintf('SVx     = %6.3f\n',SVx);
% fprintf('kappa_x = %6.3f\n',kappa__x);
% fprintf('Kx      = %6.3f\n',Bx*Cx*Dx/tyre_coeffs.FZ0);

% % Longitudinal stiffness
% Kx_vec = zeros(size(load_vec));
% for i = 1:length(load_vec)
%   [kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(0, 0, 0, load_vec(i), tyre_data);
%   Kx_vec(i) = Bx*Cx*Dx/tyre_data.Fz0;
% end
% 
% figure('Name','Kx vs Fz')
% plot(load_vec,Kx_vec,'o-')



%% Longitudinal
%% Select tyre dataset
%dataset path
data_set_path = 'dataset/';
% dataset selection and loading

%data_set = 'Hoosier_B1464run23'; % pure lateral forces
data_set = 'Hoosier_B1464run30';  % braking/traction (pure log. force) + combined

% tyre geometric data:
% Hoosier	18.0x6.0-10
% 18 diameter in inches
% 6.0 section width in inches
% tread width in inches
diameter = 18*2.56; %
Fz0 = 220;   % [N] nominal load is given
R0  = diameter/2/100; % [m] get from nominal load R0 (m) *** TO BE CHANGED ***


fprintf('Loading dataset ...')
switch data_set
    case 'Hoosier_B1464run23'
  load ([data_set_path, 'Hoosier_B1464run23.mat']); % pure lateral
  cut_start = 27760;
  cut_end   = 54500;
    case 'Hoosier_B1464run30'
  load ([data_set_path, 'Hoosier_B1464run30.mat']); % pure longitudinal
  cut_start = 19028;
  cut_end   = 37643;
  otherwise 
  error('Not found dataset: `%s`\n', data_set) ;
  
end

% select dataset portion
smpl_range = cut_start:cut_end;

fprintf('completed!\n')
%% Plot raw data

figure
tiledlayout(6,1)

ax_list(1) = nexttile; y_range = [min(min(-FZ),0) round(max(-FZ)*1.1)];
plot(-FZ)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')

ax_list(2) = nexttile; y_range = [min(min(IA),0) round(max(IA)*1.1)];
plot(IA)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(3) = nexttile; y_range = [min(min(SA),0) round(max(SA)*1.1)];
plot(SA)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Side slip')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(4) = nexttile; y_range = [min(min(SL),0) round(max(SL)*1.1)];
plot(SL)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Longitudinal slip')
xlabel('Samples [-]')
ylabel('[-]')

ax_list(5) = nexttile; y_range = [min(min(P),0) round(max(P)*1.1)];
plot(P)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Tyre pressure')
xlabel('Samples [-]')
ylabel('[psi]')

ax_list(6) = nexttile;  y_range = [min(min(TSTC),0) round(max(TSTC)*1.1)];
plot(TSTC,'DisplayName','Center')
hold on
plot(TSTI,'DisplayName','Internal')
plot(TSTO,'DisplayName','Outboard')
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Tyre temperatures')
xlabel('Samples [-]')
ylabel('[degC]')

linkaxes(ax_list,'x')


%plot(SA,FY)


%% Select some specific data
% Cut crappy data and select only 12 psi data

vec_samples = 1:1:length(smpl_range);
tyre_data = table(); % create empty table
% store raw data in table
tyre_data.SL =  SL(smpl_range);
tyre_data.SA =  SA(smpl_range)*to_rad;
tyre_data.FZ = -FZ(smpl_range);  % 0.453592  lb/kg
tyre_data.FX =  FX(smpl_range);
tyre_data.FY =  FY(smpl_range);
tyre_data.MZ =  MZ(smpl_range);
tyre_data.IA =  IA(smpl_range)*to_rad;

% Extract points at constant inclination angle
GAMMA_tol = 0.05*to_rad;
idx.GAMMA_0 = 0.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 0.0*to_rad+GAMMA_tol;
idx.GAMMA_1 = 1.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 1.0*to_rad+GAMMA_tol;
idx.GAMMA_2 = 2.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 2.0*to_rad+GAMMA_tol;
idx.GAMMA_3 = 3.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 3.0*to_rad+GAMMA_tol;
idx.GAMMA_4 = 4.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 4.0*to_rad+GAMMA_tol;
idx.GAMMA_5 = 5.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 5.0*to_rad+GAMMA_tol;
GAMMA_0  = tyre_data( idx.GAMMA_0, : );
GAMMA_1  = tyre_data( idx.GAMMA_1, : );
GAMMA_2  = tyre_data( idx.GAMMA_2, : );
GAMMA_3  = tyre_data( idx.GAMMA_3, : );
GAMMA_4  = tyre_data( idx.GAMMA_4, : );
GAMMA_5  = tyre_data( idx.GAMMA_5, : );

% Extract points at constant vertical load
% Test data done at: 
%  - 50lbf  ( 50*0.453592*9.81 =  223N )
%  - 150lbf (150*0.453592*9.81 =  667N )
%  - 200lbf (200*0.453592*9.81 =  890N )
%  - 250lbf (250*0.453592*9.81 = 1120N )

FZ_tol = 100;
idx.FZ_220  = 220-FZ_tol < tyre_data.FZ & tyre_data.FZ < 220+FZ_tol;
idx.FZ_440  = 440-FZ_tol < tyre_data.FZ & tyre_data.FZ < 440+FZ_tol;
idx.FZ_700  = 700-FZ_tol < tyre_data.FZ & tyre_data.FZ < 700+FZ_tol;
idx.FZ_900  = 900-FZ_tol < tyre_data.FZ & tyre_data.FZ < 900+FZ_tol;
idx.FZ_1120 = 1120-FZ_tol < tyre_data.FZ & tyre_data.FZ < 1120+FZ_tol;
FZ_220  = tyre_data( idx.FZ_220, : );
FZ_440  = tyre_data( idx.FZ_440, : );
FZ_700  = tyre_data( idx.FZ_700, : );
FZ_900  = tyre_data( idx.FZ_900, : );
FZ_1120 = tyre_data( idx.FZ_1120, : );

% The slip angle is varied continuously between -4 and +12° and then
% between -12° and +4° for the pure slip case

% The slip angle is varied step wise for longitudinal slip tests
% 0° , - 3° , -6 °
SA_tol = 0.5*to_rad;
idx.SA_0    =  0-SA_tol          < tyre_data.SA & tyre_data.SA < 0+SA_tol;
idx.SA_3neg = -(3*to_rad+SA_tol) < tyre_data.SA & tyre_data.SA < -3*to_rad+SA_tol;
idx.SA_6neg = -(6*to_rad+SA_tol) < tyre_data.SA & tyre_data.SA < -6*to_rad+SA_tol;
SA_0     = tyre_data( idx.SA_0, : );
SA_3neg  = tyre_data( idx.SA_3neg, : );
SA_6neg  = tyre_data( idx.SA_6neg, : );


figure()
tiledlayout(3,1)

ax_list(1) = nexttile;
plot(tyre_data.IA*to_deg)
hold on
plot(vec_samples(idx.GAMMA_0),GAMMA_0.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_1),GAMMA_1.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_2),GAMMA_2.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_3),GAMMA_3.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_4),GAMMA_4.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_5),GAMMA_5.IA*to_deg,'.');
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(2) = nexttile;
plot(tyre_data.FZ)
hold on
plot(vec_samples(idx.FZ_220),FZ_220.FZ,'.');
plot(vec_samples(idx.FZ_440),FZ_440.FZ,'.');
plot(vec_samples(idx.FZ_700),FZ_700.FZ,'.');
plot(vec_samples(idx.FZ_900),FZ_900.FZ,'.');
plot(vec_samples(idx.FZ_1120),FZ_1120.FZ,'.');
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')


ax_list(3) = nexttile;
plot(tyre_data.SA*to_deg)
hold on
plot(vec_samples(idx.SA_0),   SA_0.SA*to_deg,'.');
plot(vec_samples(idx.SA_3neg),SA_3neg.SA*to_deg,'.');
plot(vec_samples(idx.SA_6neg),SA_6neg.SA*to_deg,'.');
title('Slide slip')
xlabel('Samples [-]')
ylabel('[rad]')



%%-------------------------------------------------------------------------
% FITTING WITH GUESS VALUES and nominal vertical load
%--------------------------------------------------------------------------
%% Intersect tables to obtain specific sub-datasets

[TData0, ~] = intersect_table_data( SA_0, GAMMA_0, FZ_220 );

%% plot_selected_data

figure('Name','Selected-data')
plot_selected_data(TData0);


%% Fitting with Fz=Fz_nom= 220N and camber=0  alpha = 0 VX= 10
% ------------------
% long slip

% Fit the coeffs {pCx1, pDx1, pEx1, pEx4, pKx1, pHx1, pVx1}
FZ0 = mean(TData0.FZ);

zeros_vec = zeros(size(TData0.SL));
ones_vec  = ones(size(TData0.SL));

%Fz=Fz_nom= 220N and camber=0  alpha = 0 and k not 0
FX0_guess = MF96_FX0_vec(TData0.SL, zeros_vec, zeros_vec, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);

% check guess 
figure('Name','FX0 guess')
plot(TData0.SL,TData0.FX,'.')
hold on
plot(TData0.SL,FX0_guess,'-')


% Plot raw data and initial guess
% figure()
% plot(TDataSub.KAPPA,TDataSub.FX,'o')
% hold on
% plot(TDataSub.KAPPA,FX0_guess,'x')

% Guess values for parameters to be optimised
%    [pCx1 pDx1 pEx1 pEx4  pHx1  pKx1  pVx1 
P0 = [  1,   2,   1,  0,   0,   1,   0]; 

% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
% 1< pCx1 < 2 
% 0< pEx1 < 1 
%    [pCx1 pDx1 pEx1 pEx4  pHx1  pKx1  pVx1 
lb = [1,   0.1,   0,   0,  -10,    0,   -10];
ub = [2,    4,   1,   1,   10,   100,  10];


KAPPA_vec = TData0.SL;
FX_vec    = TData0.FX;

% check guess
SL_vec = -0.3:0.001:0.3;
FX0_fz_nom_vec = MF96_FX0_vec(SL_vec,zeros(size(SL_vec)) , zeros(size(SL_vec)), ...
                              FZ0.*ones(size(SL_vec)),tyre_coeffs);

% 
% figure
% plot(KAPPA_vec,FX_vec,'.')
% hold on
% plot(SL_vec,FX0_fz_nom_vec,'.')
% 


% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
tyre_coeffs2 = tyre_coeffs;
% alpha = 0 (pure longitudinal)
[P_fz_nom,fval,exitflag,iteration] = fmincon(@(P)resid_pure_Fx(P,FX_vec, KAPPA_vec,0,FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Update tyre data with new optimal values                             
tyre_coeffs.pCx1 = P_fz_nom(1) ; % 1
tyre_coeffs.pDx1 = P_fz_nom(2) ;  
tyre_coeffs.pEx1 = P_fz_nom(3) ;
tyre_coeffs.pEx4 = P_fz_nom(4) ;
tyre_coeffs.pHx1 = P_fz_nom(5) ; 
tyre_coeffs.pKx1 = P_fz_nom(6) ;
tyre_coeffs.pVx1 = P_fz_nom(7) ;

FX0_fz_nom_vec = MF96_FX0_vec(SL_vec,zeros(size(SL_vec)) , zeros(size(SL_vec)), ...
                              FZ0.*ones(size(SL_vec)),tyre_coeffs);

figure('Name','Fx0(Fz0)')
plot(TData0.SL,TData0.FX,'o')
hold on
%plot(TDataSub.KAPPA,FX0_fz_nom_vec,'-')
plot(SL_vec,FX0_fz_nom_vec,'-','LineWidth',2)

xlabel('$\kappa$ [-]')
ylabel('$F_{x0}$ [N]')

%% ------------------------------------------------------------------------
% FIT COEFFICIENTS WITH VARIABLE LOAD
%--------------------------------------------------------------------------
%% Fit coefficient with variable load
% extract data with variable load
[TDataDFz, ~] = intersect_table_data( SA_0, GAMMA_0 );

% long slip

% Fit the coeffs {pCx1, pDx1, pEx1, pEx4, pKx1, pHx1, pVx1}
%FZ0 = mean(TData0.FZ);

zeros_vec = zeros(size(TDataDFz.SL));
ones_vec  = ones(size(TDataDFz.SL));

% zero camber and zero side slip alpha (pure longitudinal), k = SL not 0
FX0_guess = MF96_FX0_vec(TDataDFz.SL,zeros_vec , zeros_vec, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);

% check guess 
figure()
plot(TDataDFz.SL,TDataDFz.FX,'.')
hold on
plot(TDataDFz.SL,FX0_guess,'-')

% Plot raw data and initial guess
% figure()
% plot(TDataSub.KAPPA,TDataSub.FX,'o')
% hold on
% plot(TDataSub.KAPPA,FX0_guess,'x')

% Guess values for parameters to be optimised
%    [pDx2 pEx2 pEx3 pHx2  pKx2  pKx3  pVx2] 
P0 = [  0,   0,   0,  0,   0,   0,   0]; 


% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
% 1< pCx1 < 2 
% 0< pEx1 < 1  
lb = [];
ub = [];


KAPPA_vec = TDataDFz.SL;
FX_vec    = TDataDFz.FX;
FZ_vec    = TDataDFz.FZ;

% check guess
SL_vec = -0.3:0.001:0.3;
FX0_dfz_vec = MF96_FX0_vec(SL_vec,zeros(size(SL_vec)) , zeros(size(SL_vec)), ...
                           TDataDFz.FZ,tyre_coeffs);
% 
% figure
% plot(KAPPA_vec,FX_vec,'.')
% hold on
% plot(SL_vec,FX0_dfz_vec,'.')


% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 

[P_dfz,fval,exitflag] = fmincon(@(P)resid_pure_Fx_varFz(P,FX_vec, KAPPA_vec,0,FZ_vec, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

disp(exitflag)
% Change tyre data with new optimal values                             
tyre_coeffs.pDx2 = P_dfz(1) ; % 1
tyre_coeffs.pEx2 = P_dfz(2) ;  
tyre_coeffs.pEx3 = P_dfz(3) ;
tyre_coeffs.pHx2 = P_dfz(4) ;
tyre_coeffs.pKx2 = P_dfz(5) ; 
tyre_coeffs.pKx3 = P_dfz(6) ;
tyre_coeffs.pVx2 = P_dfz(7) ;


res_FX0_dfz_vec = resid_pure_Fx_varFz(P_dfz,FX_vec,SL_vec,0 , FZ_vec,tyre_coeffs);

tmp_zeros = zeros(size(SL_vec));
tmp_ones = ones(size(SL_vec));


FX0_fz_var_vec1 = MF96_FX0_vec(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_220.FZ)*tmp_ones,tyre_coeffs);
FX0_fz_var_vec2 = MF96_FX0_vec(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_700.FZ)*tmp_ones,tyre_coeffs);
FX0_fz_var_vec3 = MF96_FX0_vec(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_900.FZ)*tmp_ones,tyre_coeffs);
FX0_fz_var_vec4 = MF96_FX0_vec(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_1120.FZ)*tmp_ones,tyre_coeffs);


figure('Name','Fx0(Fz0)')
plot(TDataDFz.SL,TDataDFz.FX,'o')
hold on
%plot(TDataSub.KAPPA,FX0_fz_nom_vec,'-')
%plot(SL_vec,FX0_dfz_vec,'-','LineWidth',2)
plot(SL_vec,FX0_fz_var_vec1,'-','LineWidth',2)
plot(SL_vec,FX0_fz_var_vec2,'-','LineWidth',2)
plot(SL_vec,FX0_fz_var_vec3,'-','LineWidth',2)
plot(SL_vec,FX0_fz_var_vec4,'-','LineWidth',2)

xlabel('$\kappa$ [-]')
ylabel('$F_{x0}$ [N]')


[kappa__x, Bx, Cx, Dx, Ex, SVx] =MF96_FX0_coeffs(0, 0, 0, mean(FZ_220.FZ), tyre_coeffs);
Calfa_vec1_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);
[kappa__x, Bx, Cx, Dx, Ex, SVx] =MF96_FX0_coeffs(0, 0, 0, mean(FZ_700.FZ), tyre_coeffs);
Calfa_vec2_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);
[kappa__x, Bx, Cx, Dx, Ex, SVx] =MF96_FX0_coeffs(0, 0, 0, mean(FZ_900.FZ), tyre_coeffs);
Calfa_vec3_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);
[kappa__x, Bx, Cx, Dx, Ex, SVx] =MF96_FX0_coeffs(0, 0, 0, mean(FZ_1120.FZ), tyre_coeffs);
Calfa_vec4_0 = magic_formula_stiffness(kappa__x, Bx, Cx, Dx, Ex, SVx);

Calfa_vec1 = MF96_CorneringStiffness(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_220.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec2 = MF96_CorneringStiffness(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_700.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec3 = MF96_CorneringStiffness(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_900.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec4 = MF96_CorneringStiffness(SL_vec,tmp_zeros ,tmp_zeros, mean(FZ_1120.FZ)*tmp_ones,tyre_coeffs);

figure('Name','C_alpha')
subplot(2,1,1)
hold on
%plot(TDataSub.KAPPA,FX0_fz_nom_vec,'-')
plot(mean(FZ_220.FZ),Calfa_vec1_0,'+','LineWidth',2)
plot(mean(FZ_700.FZ),Calfa_vec3_0,'+','LineWidth',2)
plot(mean(FZ_900.FZ),Calfa_vec4_0,'+','LineWidth',2)
plot(mean(FZ_1120.FZ),Calfa_vec2_0,'+','LineWidth',2)
legend({'$Fz_{220}$','$Fz_{700}$','$Fz_{900}$','$Fz_{1120}$'})

subplot(2,1,2)
hold on
%plot(TDataSub.KAPPA,FX0_fz_nom_vec,'-')
plot(SL_vec,Calfa_vec1,'-','LineWidth',2)
plot(SL_vec,Calfa_vec2,'-','LineWidth',2)
plot(SL_vec,Calfa_vec3,'-','LineWidth',2)
plot(SL_vec,Calfa_vec4,'-','LineWidth',2)
legend({'$Fz_{220}$','$Fz_{700}$','$Fz_{900}$','$Fz_{1120}$'})

%% ------------------------------------------------------------------------
% FIT COEFFICIENTS WITH VARIABLE CAMBER (pure longitudinal so alpha = 0, SL not 0 and gamma not 0)
%--------------------------------------------------------------------------

% extract data with variable load
[TDataGamma, ~] = intersect_table_data( SA_0, FZ_220 );

% Fit the coeffs { pDx3}

% Guess values for parameters to be optimised
P0 = [0]; 

% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
% 1< pCx1 < 2 
% 0< pEx1 < 1 
%lb = [0, 0,  0, 0,  0,  0,  0];
%ub = [2, 1e6,1, 1,1e1,1e2,1e2];
lb = [];
ub = [];


zeros_vec = zeros(size(TDataGamma.SL));
ones_vec  = ones(size(TDataGamma.SL));

KAPPA_vec = TDataGamma.SL;
GAMMA_vec = TDataGamma.IA; 
FX_vec    = TDataGamma.FX;
FZ_vec    = TDataGamma.FZ;

figure()
plot(KAPPA_vec,FX_vec);


% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
% we choose 220 [N] as load
[P_varGamma,fval,exitflag] = fmincon(@(P)resid_pure_Fx_varGamma(P,FX_vec, KAPPA_vec,GAMMA_vec,tyre_coeffs.FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Change tyre data with new optimal values                             
tyre_coeffs.pDx3 = P_varGamma(1) ; % 1

FX0_varGamma_vec = MF96_FX0_vec(KAPPA_vec,zeros_vec , GAMMA_vec, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);

figure('Name','Fx0 vs Gamma')
plot(KAPPA_vec,TDataGamma.FX,'o')
hold on
plot(KAPPA_vec,FX0_varGamma_vec,'-')
xlabel('$\kappa$ [-]')
ylabel('$F_{x0}$ [N]')
% Calculate the residuals with the optimal solution found above
res_Fx0_varGamma  = resid_pure_Fx_varGamma(P_varGamma,FX_vec, KAPPA_vec,GAMMA_vec,tyre_coeffs.FZ0, tyre_coeffs);

% R-squared is 
% 1-SSE/SST
% SSE/SST = res_Fx0_nom

% SSE is the sum of squared error,  SST is the sum of squared total
fprintf('R-squared = %6.3f\n',1-res_Fx0_varGamma);


[kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(0, 0, GAMMA_vec(3), tyre_coeffs.FZ0, tyre_coeffs);
% 
fprintf('Bx      = %6.3f\n',Bx);
fprintf('Cx      = %6.3f\n',Cx);
fprintf('mux      = %6.3f\n',Dx/tyre_coeffs.FZ0);
fprintf('Ex      = %6.3f\n',Ex);
fprintf('SVx     = %6.3f\n',SVx);
fprintf('kappa_x = %6.3f\n',kappa__x);
fprintf('Kx      = %6.3f\n',Bx*Cx*Dx/tyre_coeffs.FZ0);

% % Longitudinal stiffness
% Kx_vec = zeros(size(load_vec));
% for i = 1:length(load_vec)
%   [kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(0, 0, 0, load_vec(i), tyre_data);
%   Kx_vec(i) = Bx*Cx*Dx/tyre_data.Fz0;
% end
% 
% figure('Name','Kx vs Fz')
% plot(load_vec,Kx_vec,'o-')

%% Combined
%% Longitudinal fitting
% [TData2,~] = intersect_table_data(FZ_220, GAMMA_0, SA_3neg);
% zeros_vec = zeros(size(TData2.SL));
% ones_vec  = ones(size(TData2.SL));
% 
% KAPPA_vec = TDataGamma.SL;
% ALPHA_vec = TDataGamma.SA; 
% FX_vec    = TDataGamma.FX;
% FZ_vec    = TDataGamma.FZ;
% 
% P0 = [0.1,0.1,0.1,0.1];
% 
% [P_comb_x,fval,exitflag] = fmincon(@(P)resid_comb_Fx(P,FX_vec, KAPPA_vec, ALPHA_vec, zeros_vec,tyre_coeffs.FZ0, tyre_coeffs),...
%                                P0,[],[],[],[],lb,ub);


%% Self-Aligning
%% Select tyre dataset
%dataset path
data_set_path = 'dataset/';
% dataset selection and loading

data_set = 'Hoosier_B1464run23'; % pure lateral forces
%data_set = 'Hoosier_B1464run30';  % braking/traction (pure log. force) + combined

% tyre geometric data:
% Hoosier	18.0x6.0-10
% 18 diameter in inches
% 6.0 section width in inches
% tread width in inches
diameter = 18*2.56; %
Fz0 = 220;   % [N] nominal load is given
R0  = diameter/2/100; % [m] get from nominal load R0 (m) *** TO BE CHANGED ***


fprintf('Loading dataset ...')
switch data_set
    case 'Hoosier_B1464run23'
  load ([data_set_path, 'Hoosier_B1464run23.mat']); % pure lateral
  cut_start = 31300;
  cut_end   = 54500;
    case 'Hoosier_B1464run30'
  load ([data_set_path, 'Hoosier_B1464run30.mat']); % pure longitudinal
  cut_start = 19028;
  cut_end   = 37643;
  otherwise 
  error('Not found dataset: `%s`\n', data_set) ;
  
end

% select dataset portion
smpl_range = cut_start:cut_end;

fprintf('completed!\n')
%% Plot raw data

figure
tiledlayout(6,1)

ax_list(1) = nexttile; y_range = [min(min(-FZ),0) round(max(-FZ)*1.1)];
plot(-FZ)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')

ax_list(2) = nexttile; y_range = [min(min(IA),0) round(max(IA)*1.1)];
plot(IA)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(3) = nexttile; y_range = [min(min(SA),0) round(max(SA)*1.1)];
plot(SA)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Side slip')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(4) = nexttile; y_range = [min(min(SL),0) round(max(SL)*1.1)];
plot(SL)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Longitudinal slip')
xlabel('Samples [-]')
ylabel('[-]')

ax_list(5) = nexttile; y_range = [min(min(P),0) round(max(P)*1.1)];
plot(P)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Tyre pressure')
xlabel('Samples [-]')
ylabel('[psi]')

ax_list(6) = nexttile;  y_range = [min(min(TSTC),0) round(max(TSTC)*1.1)];
plot(TSTC,'DisplayName','Center')
hold on
plot(TSTI,'DisplayName','Internal')
plot(TSTO,'DisplayName','Outboard')
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Tyre temperatures')
xlabel('Samples [-]')
ylabel('[degC]')

linkaxes(ax_list,'x')


%plot(SA,FY)

%% Select some specific data
% Cut crappy data and select only 12 psi data
clear idx;
vec_samples = 1:1:length(smpl_range);
tyre_data = table(); % create empty table
% store raw data in table
tyre_data.SL =  SL(smpl_range);
tyre_data.SA =  SA(smpl_range)*to_rad;
tyre_data.FZ = -FZ(smpl_range);  % 0.453592  lb/kg
tyre_data.FX =  FX(smpl_range);
tyre_data.FY =  FY(smpl_range);
tyre_data.MZ =  MZ(smpl_range);
tyre_data.IA =  IA(smpl_range)*to_rad;

% Extract points at constant inclination angle
GAMMA_tol = 0.05*to_rad;
idx.GAMMA_0 = 0.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 0.0*to_rad+GAMMA_tol;
idx.GAMMA_1 = 1.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 1.0*to_rad+GAMMA_tol;
idx.GAMMA_2 = 2.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 2.0*to_rad+GAMMA_tol;
idx.GAMMA_3 = 3.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 3.0*to_rad+GAMMA_tol;
idx.GAMMA_4 = 4.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 4.0*to_rad+GAMMA_tol;
idx.GAMMA_5 = 5.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 5.0*to_rad+GAMMA_tol;
GAMMA_0  = tyre_data( idx.GAMMA_0, : );
GAMMA_1  = tyre_data( idx.GAMMA_1, : );
GAMMA_2  = tyre_data( idx.GAMMA_2, : );
GAMMA_3  = tyre_data( idx.GAMMA_3, : );
GAMMA_4  = tyre_data( idx.GAMMA_4, : );
GAMMA_5  = tyre_data( idx.GAMMA_5, : );

% Extract points at constant vertical load
% Test data done at: 
%  - 50lbf  ( 50*0.453592*9.81 =  223N )
%  - 150lbf (150*0.453592*9.81 =  667N )
%  - 200lbf (200*0.453592*9.81 =  890N )
%  - 250lbf (250*0.453592*9.81 = 1120N )

FZ_tol = 100;
idx.FZ_220  = 220-FZ_tol < tyre_data.FZ & tyre_data.FZ < 220+FZ_tol;
idx.FZ_440  = 440-FZ_tol < tyre_data.FZ & tyre_data.FZ < 440+FZ_tol;
idx.FZ_700  = 700-FZ_tol < tyre_data.FZ & tyre_data.FZ < 700+FZ_tol;
idx.FZ_900  = 900-FZ_tol < tyre_data.FZ & tyre_data.FZ < 900+FZ_tol;
idx.FZ_1120 = 1120-FZ_tol < tyre_data.FZ & tyre_data.FZ < 1120+FZ_tol;
FZ_220  = tyre_data( idx.FZ_220, : );
FZ_440  = tyre_data( idx.FZ_440, : );
FZ_700  = tyre_data( idx.FZ_700, : );
FZ_900  = tyre_data( idx.FZ_900, : );
FZ_1120 = tyre_data( idx.FZ_1120, : );


figure()
tiledlayout(3,1)

ax_list(1) = nexttile;
plot(tyre_data.IA*to_deg)
hold on
plot(vec_samples(idx.GAMMA_0),GAMMA_0.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_1),GAMMA_1.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_2),GAMMA_2.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_3),GAMMA_3.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_4),GAMMA_4.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_5),GAMMA_5.IA*to_deg,'.');
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(2) = nexttile;
plot(tyre_data.FZ)
hold on
plot(vec_samples(idx.FZ_220),FZ_220.FZ,'.');
plot(vec_samples(idx.FZ_440),FZ_440.FZ,'.');
plot(vec_samples(idx.FZ_700),FZ_700.FZ,'.');
plot(vec_samples(idx.FZ_900),FZ_900.FZ,'.');
plot(vec_samples(idx.FZ_1120),FZ_1120.FZ,'.');
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')


ax_list(3) = nexttile;
plot(tyre_data.SA)
hold on
title('Slide slip')
xlabel('Samples [-]')
ylabel('[rad]')

%% Intersect tables to obtain specific sub-datasets and plot them
[TData0, ~] = intersect_table_data(GAMMA_0, FZ_220 );

%%-------------------------------------------------------------------------
% FITTING WITH GUESS VALUES and nominal vertical load
%--------------------------------------------------------------------------
  tyre_coeffs = initialise_tyre_data(R0, Fz0);

%% MISSING scaling_factors: see maple MF96_MZ0_coeffs_eqns
tyre_coeffs.LKY = 1;
tyre_coeffs.LT = 1;
tyre_coeffs.LMR = 1;

%%
FY_vec = TData0.FY;
FZ0 = mean(TData0.FZ);
zeros_vec = zeros(size(TData0.SA));
ones_vec  = ones(size(TData0.SA));

% based on pure lateral Fy0(k = 0, alpha, gamma = 0)
FY0_guess = MF96_FY0_vec(zeros_vec, TData0.SA, zeros_vec, tyre_coeffs.FZ0 * ones_vec, tyre_coeffs);
% test with combined lateral Fy(k, alpha, gamma = 0 (1st case))
%FY_guess = MF96_FY_vec(zeros_vec, TData0.SA, zeros_vec, tyre_coeffs.FZ0 * ones_vec, tyre_coeffs);

MZ0_guess = MF96_MZ0_vec(zeros_vec, TData0.SA, zeros_vec, tyre_coeffs.FZ0 * ones_vec, FY0_guess, tyre_coeffs); 

% figure 4
% check guess 
figure('Name','MZ0 guess')
plot(TData0.SA, TData0.MZ,'.')
hold on
plot(TData0.SA, MZ0_guess, '-')

%% ------------------------------------------------------------------------
% FITTING WITH NOMINAL LOAD Fz=Fz_nom= 220N and camber=0  alpha = 0 VX= 10
%--------------------------------------------------------------------------
% Guess values for parameters to be optimised
%    [qHz1,qBz1,qCz1,qDz1,qEz1,qEz4,qBz9,qDz6,qBz10]
P0 = [0,    0,    0,   0,  0,   0,   0,   0,   0]; 

% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
%    [qHz1,qBz1,qCz1,qDz1,qEz1,qEz4,qBz9,qDz6,qBz10]
lb = [0, -5, -5, 0, -5, 0, 0, 0, 0, 0, -5]; % lower bound
ub = [15, 10, 10 , 10, 10, 10, 10, 50, 10 ]; % upper bound
% lb = [];
% ub = [];

ALPHA_vec = TData0.SA;  % lateral slip angle
MZ_vec    = TData0.MZ;  % aligning moment
FY_vec    = TData0.FY;  % pure lateral force

% check guess
SA_vec = -0.3:0.001:0.3; % lateral slip to be used in the plot for check the interpolation outside the input data range

% resid_pure_Mz returns the residual, so minimize the residual varying Y. 
% It is an unconstrained minimization problem 
[P_fz_nom, fval, exitflag] = fmincon(@(P)resid_pure_Mz(P, MZ_vec, ALPHA_vec, 0, FZ0, FY_vec, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

P_fz_nom
fval

% Update tyre data with new optimal values                   
tyre_coeffs.qHz1  = P_fz_nom(1); 
tyre_coeffs.qBz1  = P_fz_nom(2);
tyre_coeffs.qCz1  = P_fz_nom(3);
tyre_coeffs.qDz1  = P_fz_nom(4);
tyre_coeffs.qEz1  = P_fz_nom(5);
tyre_coeffs.qEz4  = P_fz_nom(6);
tyre_coeffs.qBz9  = P_fz_nom(7);
tyre_coeffs.qDz6  = P_fz_nom(8);
tyre_coeffs.qBz10 = P_fz_nom(9);

FY0_fz_nom_vec = MF96_FY0_vec(zeros(size(SA_vec)), SA_vec, zeros(size(SA_vec)), tyre_coeffs.FZ0.*ones(size(SA_vec)), tyre_coeffs);
MZ0_fz_nom_vec = MF96_MZ0_vec(zeros(size(SA_vec)), SA_vec, zeros(size(SA_vec)), FZ0.*ones(size(SA_vec)), FY0_fz_nom_vec, tyre_coeffs);
        
% figure 5          
plot(TData0.SA, TData0.MZ,'o')
hold on
plot(SA_vec, MZ0_fz_nom_vec, '-', 'LineWidth', 2)
xlabel('$\alpha$ [rad]')
ylabel('$M_{z0}$ [N m]')

%% ------------------------------------------------------------------------
% FIT COEFFICIENTS WITH VARIABLE LOAD
%--------------------------------------------------------------------------
% extract data with variable load
[TDataDFz, ~] = GAMMA_0; % intersect_table_data(SA_0, GAMMA_0);

zeros_vec = zeros(size(TDataDFz.SA));
ones_vec  = ones(size(TDataDFz.SA));

FY0_guess = MF96_FY0_vec(zeros_vec, GAMMA_0.SA, zeros_vec, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);

MZ0_guess = MF96_MZ0_vec(zeros_vec, GAMMA_0.SA, zeros_vec, tyre_coeffs.FZ0*ones_vec, FY0_guess, tyre_coeffs);

% check guess 
figure()
plot(GAMMA_0.SA,GAMMA_0.MZ,'.')
hold on
plot(GAMMA_0.SA,MZ0_guess,'-')

% Guess values for parameters to be optimised
%    [  qHz2 qHz3 qHz4 qBz2 qBz3 qEz2 qEz3 qDz7 qDz8 qDz9]
%% WORK HERE
% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
% 1< pCx1 < 2 
% 0< pEx1 < 1 
%    [ qHz2, qBz2, qBz3, qDz2, qEz2, qEz3, qDz7]
P0 = [0, 0, 0, 0, 0, 0, 0];
%% WORK HERE
lb = [0, 0, 0, 0, 0, 0, 0];
ub = [10, 10, 10, 10, 10, 10, 10];

ALPHA_vec = TDataDFz.SA;
MZ_vec    = TDataDFz.MZ;
FZ_vec    = TDataDFz.FZ;
FY_vec    = TDataDFz.FY;

% check guess
SA_vec = -0.3:0.001:0.3;
% FX0_dfz_vec = MF96_FX0_vec(SL_vec,zeros(size(SL_vec)), zeros(size(SL_vec)), ...
%                            TDataDFz.FZ,tyre_coeffs);

% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
[P_dfz,fval,exitflag] = fmincon(@(P)resid_pure_Mz_varFz(P, MZ_vec, ALPHA_vec, 0, FZ_vec, FY_vec, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

disp(exitflag)
% Change tyre data with new optimal values          
tyre_coeffs.qHz2 = P(1); 
tyre_coeffs.qBz2 = P(2); 
tyre_coeffs.qBz3 = P(3); 
tyre_coeffs.qDz2 = P(4); 
tyre_coeffs.qEz2 = P(5); 
tyre_coeffs.qEz3 = P(6); 
tyre_coeffs.qDz7 = P(7); 


res_MZ0_dfz_vec = resid_pure_Mz_varFz(P_dfz,MZ_vec,SA_vec, 0, FZ_vec,FY_vec,tyre_coeffs);

tmp_zeros = zeros(size(SA_vec));
tmp_ones = ones(size(SA_vec));
%%

FY0_fz_var_vec1 =  MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_220.FZ)*tmp_ones, tyre_coeffs);
MZ0_fz_var_vec1 = MF96_MZ0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_220.FZ)*tmp_ones, FY0_fz_var_vec1, tyre_coeffs);

FY0_fz_var_vec2 =  MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_440.FZ)*tmp_ones, tyre_coeffs);
MZ0_fz_var_vec2 = MF96_MZ0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_440.FZ)*tmp_ones, FY0_fz_var_vec2, tyre_coeffs);

FY0_fz_var_vec3 =  MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_900.FZ)*tmp_ones, tyre_coeffs);
MZ0_fz_var_vec3 = MF96_MZ0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_900.FZ)*tmp_ones, FY0_fz_var_vec3, tyre_coeffs);

FY0_fz_var_vec4 =  MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_1120.FZ)*tmp_ones, tyre_coeffs);
MZ0_fz_var_vec4 = MF96_MZ0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_1120.FZ)*tmp_ones, FY0_fz_var_vec4, tyre_coeffs);

FY0_fz_var_vec5 =  MF96_FY0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_700.FZ)*tmp_ones, tyre_coeffs);
MZ0_fz_var_vec5 = MF96_MZ0_vec(tmp_zeros, SA_vec, tmp_zeros, mean(FZ_700.FZ)*tmp_ones, FY0_fz_var_vec5, tyre_coeffs);


figure('Name','Mz0(Fz0)')
plot(GAMMA_0.SA,GAMMA_0.MZ,'o')
hold on
plot(SA_vec,MZ0_fz_var_vec1,'-','LineWidth',2)
plot(SA_vec,MZ0_fz_var_vec2,'-','LineWidth',2)
plot(SA_vec,MZ0_fz_var_vec3,'-','LineWidth',2)
plot(SA_vec,MZ0_fz_var_vec4,'-','LineWidth',2)
plot(SA_vec,MZ0_fz_var_vec5,'-','LineWidth',2)

xlabel('$\alpha$ [rad]')
ylabel('$M_{z0}$ [N m]')


%% Stiffness
%[alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(0, 0, mean(FZ_220.FZ), tyre_coeffs);
% Calfa_vec1_0 = magic_formula_stiffness_MZ(alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr);
%[alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(0, 0, mean(FZ_700.FZ), tyre_coeffs);
% Calfa_vec2_0 = magic_formula_stiffness_MZ(alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr);
%[alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(0, 0, mean(FZ_900.FZ), tyre_coeffs);
% Calfa_vec3_0 = magic_formula_stiffness_MZ(alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr);
%[alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(0, 0, mean(FZ_1120.FZ), tyre_coeffs);
% Calfa_vec4_0 = magic_formula_stiffness_MZ(alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr);

%Calfa_vec1 = MF96_CorneringStiffness_MZ(SA_vec, tmp_zeros, mean(FZ_220.FZ)  * tmp_ones, FY0_fz_var_vec1, tyre_coeffs);
%Calfa_vec1 = MF96_CorneringStiffness_MZ(SA_vec, tmp_zeros, mean(FZ_440.FZ)  * tmp_ones, FY0_fz_var_vec2, tyre_coeffs);
%Calfa_vec3 = MF96_CorneringStiffness_MZ(SA_vec, tmp_zeros, mean(FZ_700.FZ)  * tmp_ones, FY0_fz_var_vec3, tyre_coeffs);
%Calfa_vec4 = MF96_CorneringStiffness_MZ(SA_vec, tmp_zeros, mean(FZ_1120.FZ) * tmp_ones, FY0_fz_var_vec4, tyre_coeffs);

% figure 7
% plot_stiffness_MZ; 

%% ------------------------------------------------------------------------
% FIT COEFFICIENTS WITH VARIABLE CAMBER
%--------------------------------------------------------------------------
% extract data with variable load
[TDataGamma, ~] = FZ_220; % intersect_table_data(SA_0, FZ_220);

% Guess values for parameters to be optimised
%    [qHz3, qBz4, qBz5, qDz3, qDz4, qEz5, qDz8]
P0 = [1, 1, 1, 1, 1, 1, 1]; 

% NOTE: many local minima => limits on parameters are fundamentals
lb = [0, 0, 0, 0, 0, 0, 0];
ub = [10, 10, 10, 10, 10, 10, 10];

zeros_vec = zeros(size(TDataGamma.SA));
ones_vec  = ones(size(TDataGamma.SA));

ALPHA_vec = TDataGamma.SA;
GAMMA_vec = TDataGamma.IA; 
MZ_vec    = TDataGamma.MZ;
FZ_vec    = TDataGamma.FZ;
FY_vec    = TDataGamma.FY;

figure()
plot(ALPHA_vec,MZ_vec);

TDataGamma0 = intersect_table_data(GAMMA_0, FZ_220);
TDataGamma1 = intersect_table_data(GAMMA_1, FZ_220);
TDataGamma2 = intersect_table_data(GAMMA_2, FZ_220);
TDataGamma3 = intersect_table_data(GAMMA_3, FZ_220);
TDataGamma4 = intersect_table_data(GAMMA_4, FZ_220);

FY0_Gamma0 = MF96_FY0_vec(zeros_vec, TDataGamma0.SA, TDataGamma0.IA, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);
MZ0_Gamma0 = MF96_MZ0_vec(TDataGamma0.SA, TDataGamma0.IA, tyre_coeffs.FZ0*ones_vec, FY0_Gamma0, tyre_coeffs);

FY0_Gamma1 = MF96_FY0_vec(zeros_vec, TDataGamma1.SA, TDataGamma1.IA, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);
MZ0_Gamma1 = MF96_MZ0_vec(zeros_vec, TDataGamma1.SA, TDataGamma1.IA, tyre_coeffs.FZ0*ones_vec, FY0_Gamma1, tyre_coeffs);

FY0_Gamma2 = MF96_FY0_vec(zeros_vec, TDataGamma2.SA, TDataGamma2.IA, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);
MZ0_Gamma2 = MF96_MZ0_vec(zeros_vec, TDataGamma2.SA, TDataGamma2.IA, tyre_coeffs.FZ0*ones_vec, FY0_Gamma2, tyre_coeffs);

FY0_Gamma3 = MF96_FY0_vec(zeros_vec, TDataGamma3.SA, TDataGamma3.IA, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);
MZ0_Gamma3 = MF96_MZ0_vec(zeros_vec, TDataGamma3.SA, TDataGamma3.IA, tyre_coeffs.FZ0*ones_vec, FY0_Gamma3, tyre_coeffs);

FY0_Gamma4 = MF96_FY0_vec(zeros_vec, TDataGamma4.SA, TDataGamma4.IA, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);
MZ0_Gamma4 = MF96_MZ0_vec(zeros_vec, TDataGamma4.SA, TDataGamma4.IA, tyre_coeffs.FZ0*ones_vec, FY0_Gamma4, tyre_coeffs);

plot(TDataGamma0.SA, TDataGamma0.MZ)
hold on
plot(TDataGamma0.SA, MZ0_Gamma0)

plot(TDataGamma1.SA, TDataGamma1.MZ);
hold on
plot(TDataGamma1.SA, MZ0_Gamma1);

plot(TDataGamma2.SA, TDataGamma2.MZ);
hold on
plot(TDataGamma2.SA, MZ0_Gamma2);

plot(TDataGamma3.SA, TDataGamma3.MZ);
hold on
plot(TDataGamma3.SA, MZ0_Gamma3);

plot(TDataGamma4.SA, TDataGamma4.MZ)
hold on
plot(TDataGamma4.SA, MZ0_Gamma4)
legend('location', 'northeast')
xlabel('$\alpha$')
ylabel('$M_{Z0}$ [N]')
title('Fitting with variable camber, $F_{Z}$ = 220[N]', 'FontSize',font_size_title)
grid on
%export_fig(fig_fit_variable_camber_MZ, 'images\fig_fit_variable_camber_MZ.png');

% fig_camber_MZ = figure('Color', 'w');
% plot(ALPHA_vec, MZ_vec, '.', 'DisplayName', '$F_{z}=220 [N]$');
% grid on
% xlabel('$\gamma$ [rad]')
% ylabel('$M_Z$ [N]')
% title('Self aligning moment as function of camber angle')
% legend('Location', 'best')
% export_fig(fig_camber_MZ, 'images\fig_camber_MZ.png')

% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 
[P_varGamma, fval, exitflag] = fmincon(@(P)resid_pure_Mz_varGamma(P, MZ_vec, ALPHA_vec, GAMMA_vec, tyre_coeffs.FZ0, FY_vec, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Change tyre data with new optimal values
tyre_coeffs.qBz4 = P_varGamma(1); 
tyre_coeffs.qBz5 = P_varGamma(2); 
tyre_coeffs.qEz5 = P_varGamma(3); 
%%
%% FY0_varGamma_vec = MF96_FY0_vec(zeros_vec, ALPHA_vec, GAMMA_vec, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);
%% MZ0_varGamma_vec = MF96_MZ0_vec(zeros_vec, ALPHA_vec, GAMMA_vec, tyre_coeffs.FZ0*ones_vec, FY0_varGamma_vec, tyre_coeffs);

plot(TDataGamma.SA, TDataGamma.MZ) 
xlabel('$\alpha$ [rad]')
ylabel('$M_{z0}$ [N m]')


FY0_varGamma_vec0 = MF96_FY0_vec(zeros_vec, SA_vec, GAMMA_0.IA, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);
MZ0_varGamma_vec0 = MF96_MZ0_vec(zeros_vec, SA_vec, GAMMA_0.IA, tyre_coeffs.FZ0*ones_vec,FY0_varGamma_vec0,tyre_coeffs);

FY0_varGamma_vec1 = MF96_FY0_vec(zeros_vec, SA_vec, GAMMA_1.IA, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);
MZ0_varGamma_vec1 = MF96_MZ0_vec(zeros_vec, SA_vec, GAMMA_1.IA, tyre_coeffs.FZ0*ones_vec,FY0_varGamma_vec1,tyre_coeffs);

FY0_varGamma_vec2 = MF96_FY0_vec(zeros_vec, SA_vec, GAMMA_2.IA, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);
MZ0_varGamma_vec2 = MF96_MZ0_vec(zeros_vec, SA_vec, GAMMA_2.IA, tyre_coeffs.FZ0*ones_vec,FY0_varGamma_vec2,tyre_coeffs);

FY0_varGamma_vec3 = MF96_FY0_vec(zeros_vec, SA_vec, GAMMA_3.IA, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);
MZ0_varGamma_vec3 = MF96_MZ0_vec(zeros_vec, SA_vec, GAMMA_3.IA, tyre_coeffs.FZ0*ones_vec,FY0_varGamma_vec3,tyre_coeffs);

FY0_varGamma_vec4 = MF96_FY0_vec(zeros_vec, SA_vec, GAMMA_4.IA, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);
MZ0_varGamma_vec4 = MF96_MZ0_vec(zeros_vec, SA_vec, GAMMA_4.IA, tyre_coeffs.FZ0*ones_vec,FY0_varGamma_vec4,tyre_coeffs);

figure('Name','Fy0 vs Gamma')
plot(ALPHA_vec,FZ_220.MZ,'o')
hold on
plot(SA_vec,MZ0_varGamma_vec0,'-')
plot(SA_vec,MZ0_varGamma_vec1,'-')
plot(SA_vec,MZ0_varGamma_vec2,'-')
plot(SA_vec,MZ0_varGamma_vec3,'-')
plot(SA_vec,MZ0_varGamma_vec4,'-')
xlabel('$\alpha$ [rad]')
ylabel('$M_{z0}$ [N m]')

% R-squared is 
% 1-SSE/SST
% SSE/SST = res_Fx0_nom

% SSE is the sum of squared error,  SST is the sum of squared total
fprintf('R-squared = %6.3f\n', 1 - res_Mz0_varGamma);

% [~, ~, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(alpha, phi, Fz, tyre_data)
% fprintf('By      = %6.3f\n', By);
% fprintf('Cy      = %6.3f\n', Cy);
% fprintf('muy     = %6.3f\n', Dy/tyre_coeffs.FZ0);
% fprintf('Ey      = %6.3f\n', Ey);
% fprintf('SVy     = %6.3f\n', SVy);
% fprintf('alpha_y = %6.3f\n', alpha__y);
% fprintf('Ky      = %6.3f\n', By*Cy*Dy/tyre_coeffs.FZ0);


%% Save tyre data structure to mat file
%
save(['tyre_' data_set,'.mat'],'tyre_coeffs');


