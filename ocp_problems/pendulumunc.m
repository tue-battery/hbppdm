function [OCPval] = pendulumc(N)
%% Inverted pendulum continuous-time model
%  [x;
%   x_dot;
%   theta;
%   theta_dot]
%

Ts = 0.1;
xinit = [4;-5;0;0];

% Parameters inverted pendulum
m = 0.2;
M = 1;
b = 0.05;
I = 0.01;
g = 9.8;
l = 0.5;
p = (I+m*l^2)*(M+m)-m^2*l^2;
Ac = [0      1              0           0; 
      0 -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0; 
      0      0              0           1;  
      0 -(m*l*b)/p       m*g*l*(M+m)/p  0];
Bc = [     0, 0;    
    (I+m*l^2)/p, m*l/p;  
    0,0;   
    m*l/p,(I+m*l^2)/p];
clear m M b I g l p;
%Discretization 
% tau = 0.001;
ssd = c2d(ss(Ac,Bc,eye(4),0),Ts);
Ad = ssd.A;
Bd = ssd.B;
Cd = ssd.C;

%Cost Matrices
Q = diag([1 1 1 1]);
R = 0.1*eye(2);

%Teminal Cost
P=idare(ssd.A,ssd.B,Q,R);

% Bounds

uubt = [4;4];
ulbt = [-4;-4];

nStates = 4;
nInputs = 2;

time = (0:Ts:floor(0.9*N)*Ts)';
fs = 1/9;
%xubs = [(3.00*sin(fs*2*pi*time+0.5*pi)+1); 5*ones(N-floor(0.9*N),1)];

xub = [10*ones(N+1,1),5*ones(N+1,1),10*ones(N+1,1),5*ones(N+1,1)];
xlb = -[10*ones(N+1,1),5*ones(N+1,1),10*ones(N+1,1),5*ones(N+1,1)];
uub = repmat(uubt,1,N);
ulb = repmat(ulbt,1,N);
F = zeros(nStates,N+1);

OCPval = genOCPproblem(ssd,Q,P,R,xub,xlb,uubt,ulbt,N,xinit);
end

