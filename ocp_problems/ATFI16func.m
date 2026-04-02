function [OCPval] = ATFI16func(N)
Ts = 0.05; % Sampling time
xinit = [0;0;0;0];
% plant

A = [-0.0151 -60.5651 0 -32.1740; -0.0001 -1.3411 0.9929 0; 0.0002 43.2541 -0.8694 0; 0 0 1 0];
%B = [-0.08 -0.635; -0.029 -0.014; -0.868 -0.092; -0.022 -0.002];
B = [-2.516 -13.136 ; -.1689 -.2514 ; -17.251 -1.5766 ; 0 0];
C = [0 1 0 0; 0 0 0 1];

nStates = 4;
nInputs = 2;

ssc = ss(A,B,C,0);
ssd = c2d(ssc,Ts);

% matrices for solving ocp
Qy = 1e2*eye(2);
Qx = diag([1e-4 0 1e-3 0]);
Q = ssd.C'*Qy*ssd.C+Qx;
P = Q;
R = 2*1e-2*eye(2);

% constraints
time = (0:Ts:N*Ts)';
fs = 1/7;
xub = ones(length(time),4).*[1e6,0.5,1e6,100];

xlb = -xub;%ones(length(time),4).*-[1000,0.5,1000,100];

uub = [25;25];
ulb = -[25;25];

OCPval = genOCPproblem(ssd,Q,P,R,xub,xlb,uub,ulb,N,xinit);

r = [zeros(OCPval.N+1,1),[10*ones(ceil(length(time)/2),1);zeros(floor(length(time)/2),1)]];

for k = 1:OCPval.N+1
%    OCPval.F((k-1)*OCPval.nStates+1:k*OCPval.nStates,1) = (r(k,:)*Qy*ssd.C)';
   OCPval.F((k-1)*OCPval.nStates+1:k*OCPval.nStates,1) = -ssd.C'*Qy*r(k,:)';
end
end

