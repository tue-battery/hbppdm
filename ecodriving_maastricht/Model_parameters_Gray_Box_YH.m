function [veh_par, EM, bat, diesel, scl] = Model_parameters_Gray_Box_YH()
    % Vehicle
    veh_par.m = 14375;                  % Vehicle mass [kg]
    % veh_par.m = 1437;                 % Vehicle mass [kg]
    veh_par.Cd = 2.733827413043453;     % Total value for air resistance [-]
    veh_par.Cr = 0.006274889551216;     % vehicle friction drag [-]
    veh_par.mu = 0.35;                  % friction coefficient 
    veh_par.g = 9.81;                   % gravitational force [m/s^2]
    veh_par.a_Min = -1.5;               % maximum deceleration [m/s^2]
    veh_par.a_Max = 1.5;                % maximum acceleration [m/s^2]
    veh_par.lr = 1.003414231681762e-04; % Distance from rear axle to gravity center [m]

    % EM
    EM.C_m = [5.452147502873735e-04	-3.355660822579786	0.693450490035610;
        0	76.098384230524870	-4.224720576801110;
        0	0	5.644443339616729]; % EM fitting matrix (Cholesky's)

    EM.PdCm1 = 62.9510606673657; % Force limit coefficient 1
    EM.PdCm2 = 156199.072522853; % Force limit coefficient 2
    EM.PtCm = 33255.6000000000;  % Force limit 
    
    
    % Battery
    bat.Q = -6.156497431526837e-08; % Quadratic fit. coef. 
    bat.F = 1.000097195829144;      % Linear fit. coef. 
    bat.Pblim = 0.83*350e3;         % Power limit battery [W]
    bat.Cap = 350*3.6e6;            % Battery capacity [J]
    
    % Diesel generator
    diesel.a = 0.4078;       % Linear fit. coef. 
    diesel.b = -3.6044e3;    % Const. fit. coef.       
    diesel.Plow = 8e3;       % Minimal power [W]
    diesel.Phigh = 145e3;    % Maximal power [W]
    
    % Scaling variables
    %scl.E_kin = 1/2*veh_par.m*max(v_max);
    %scl.EkinDyn = scl.E_kin*50;
    %scl.Friction = max(veh_par.m*veh_par.g*sin(alpha) + veh_par.Cr*veh_par.m*veh_par.g*cos(alpha) + 2/veh_par.m*scl.E_kin);
    %scl.v = max(v_max);
    %scl.LTG = 1/min(v_max);
    scl.tDyn = 20;
    scl.F_M = 60e3;
    scl.F_ME = 75e3;
    scl.F_M_lim1 = 3e5;
    scl.F_M_lim2 = 5e5;
    %scl.F_M_lim3 = 1e3;
end

