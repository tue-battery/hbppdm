function [OCPval] = eco_night_attempt(OCPval, w, N)
%eco_night_attempt Linearizes the NL dynamics for use in the solver.
% States: x = Ekin
% Inputs: u = [Fm, ell, s_Ek]

Ts = 20; % [m]
xinit = 1; % 1 MJ
nStates = 1;
nInputs = 3;
nVar    = nStates + nInputs;
constraints = 2;

sigma_Ek = 0;
sigma_fp = 15;
sigma_ell= 100;

x_scale = 1e-6;

veh_par.m = 15000; %14375;                  % Vehicle mass [kg]
veh_par.Cd = 3.11; %1.225*0.577*0.5*8.7975;% Total value for air resistance [-]
veh_par.Cr = 0.0065;                % vehicle friction drag [-]
veh_par.mu = 0.35;                  % friction coefficient 
veh_par.g = 9.81;                   % gravitational force [m/s^2]
veh_par.Fmax = 150000*x_scale;      % Newtons
veh_par.vmax = 100/3.6;             % Newtons
veh_par.vmin = 5/3.6;               % Newtons

% model with (Fme - Fm)ell = [Fm ell 1] Q [Fm ell 1]T
veh_par.Qc = [5.452147502873735e-04	-3.355660822579786	0.693450490035610;
             0	                    76.098384230524870	-4.224720576801110;
             0	                    0	                5.644443339616729];
veh_par.Q  = veh_par.Qc' * veh_par.Qc;

%Scaling to mega newton:
veh_par.Q(1,1)     = veh_par.Q(1,1)*1e6;
veh_par.Q(2:3,2:3) = veh_par.Q(2:3,2:3).* 1e-6;

% Some calculations we only do once.
if isempty(OCPval)
    
    % bounds on states
    xub =  x_scale*.5*veh_par.m*veh_par.vmax.^2 ;
    xlb =  x_scale*.5*veh_par.m*veh_par.vmin.^2 ;
    % bounds on inputs
    uubt = [ veh_par.Fmax;  1./veh_par.vmin;  inf];
    ulbt = [-veh_par.Fmax;  1./veh_par.vmax;  0];
    wub = [repmat([xub;uubt],N,1); xub];
    wlb = [repmat([xlb;ulbt],N,1); xlb];
    
    OCPval.wub = wub;
    OCPval.wlb = wlb;
    
    height = 150*sin((1:Ts:N*Ts)/.0001);
    d_height = diff(height);
    alpha = [atan(d_height/Ts)';0];
    
    FV  = (veh_par.m*veh_par.g*veh_par.Cr*cos(alpha) + veh_par.m*veh_par.g*sin(alpha))*x_scale*Ts;
    
    OCPval.FV = FV;
    OCPval.xinit = xinit;
    OCPval.N = N;
    OCPval.nInputs = nInputs;
    OCPval.nStates = nStates;
    OCPval.Ts = Ts;
end

FV = OCPval.FV;

%% Building Cost function

% Linearization point.
Ekin = [w(1:nVar:end)];
Fm  = w(2:nVar:end);
ell  = w(3:nVar:end);
v = 1./ell;

% Cost function EM
for k = 1:N
Cost_Q(:,:,k) = [0, 0, 0, 0;
                 0, v(k)*veh_par.Q(1,1), v(k)*veh_par.Q(1,2), 0;
                 0, v(k)*veh_par.Q(1,2), v(k)*veh_par.Q(2,2), 0;
                 0, 0, 0, 1e2];
Cost_F(1+(k-1)*nVar:k*nVar,1) = [0;
                               (1 + 2*v(k)*veh_par.Q(1,3));
                               2*v(k)*veh_par.Q(2,3);
                               0];
end

[gamma_row,gamma_col,gamma_val] = gen_diagonal(Cost_Q,N,0,0);
G_em = sparse(gamma_row,gamma_col,gamma_val,nStates*(N+1)+nInputs*N,nStates*(N+1)+nInputs*N);
F_em = [Cost_F;0];

% Cost function regularization
G_reg = 1*speye((N+1)*nStates+N*nInputs);
F_reg = -2*w'*G_reg;


% Cost linear
Fu = [sigma_fp; sigma_ell; 0];  
F_cost = [repmat([sigma_Ek;Fu],N,1);sigma_Ek];

% Final cost function

Gfull = G_em + G_reg;

% diagonalize
Gdiag    = sparse(1:length(Gfull),1:length(Gfull), diag(Gfull));
Gnondiag = Gfull-Gdiag;

OCPval.G  = Gdiag;
OCPval.iG = sparse(1:length(Gfull),1:length(Gfull), 1./diag(Gdiag));
OCPval.F  = F_reg' + F_em + F_cost + Gnondiag*w;



%% Equality Constraints

%Discretization 
A_d  = (1 - 2*Ts*veh_par.Cd/veh_par.m); % continuous dynamics
B_d  = [Ts 0 0];
    
%Formulating the dynamics
for k = 1:N
    Cx(:,:,k) = [A_d, B_d;
                 ell(k).^2, 0, 2*Ekin(k+1)*ell(k), -1];
    b(1+(k-1)*constraints:k*constraints,1) = [
        FV(k)
        2*Ekin(k+1)*(ell(k)^2) + 0.5*veh_par.m*x_scale];
end

% All dynamics on the diagonal
[gamma_row,gamma_col,gamma_val] = gen_diagonal(Cx,N,1,0);

% initial condition and identity for Ekin dyn
gamma_row = [gamma_row; 1; (2:constraints:(constraints*N+1))'];
gamma_col = [gamma_col; 1; (nVar+1:nVar:(nVar*N+1))']; 
gamma_val = [gamma_val; 1; -ones(N,1)];




OCPval.AA = sparse([gamma_row],[gamma_col],[gamma_val]);
OCPval.bb = [xinit;b];

% rescaling constraints.
constraint_scale = 1;%normest(OCPval.G) / normest(OCPval.AA);
OCPval.AA = OCPval.AA * constraint_scale;
OCPval.bb = OCPval.bb * constraint_scale;

end

function [gamma_row,gamma_col,gamma_val] = gen_diagonal(A,N,row_offset_A,col_offset_A)
    
    counter = 0;
    A_size = sum(A~=0,"all");

    gamma_row = zeros(A_size,1);
    gamma_col = zeros(A_size,1);
    gamma_val = zeros(A_size,1);

    for k = 1:N
        [A_row,A_col,A_var] = find(A(:,:,k));
        
        [len_A_row, len_A_col, ~] = size(A(:,:,k));
        
        len_A = length(A_row);
        
        
        gamma_row(1+counter:counter+len_A,1) = A_row +(k-1)*len_A_row +row_offset_A;
        gamma_col(1+counter:counter+len_A,1) = A_col +(k-1)*len_A_col +col_offset_A;
        gamma_val(1+counter:counter+len_A,1) = A_var;
        
        counter = counter + len_A;
    end

end