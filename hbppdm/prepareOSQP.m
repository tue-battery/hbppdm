function [m] = prepareOSQP(OCPval)
% build required matrices
A = [OCPval.AA; eye((OCPval.N+1)*OCPval.nStates), zeros((OCPval.N+1)*OCPval.nStates, OCPval.N*OCPval.nInputs); zeros(OCPval.N*OCPval.nInputs, (OCPval.N+1)*OCPval.nStates), eye(OCPval.N*OCPval.nInputs) ];
l = [OCPval.bb; OCPval.wlb];
u = [OCPval.bb; OCPval.wub];




% formulate osqp problem
m=osqp;

m.setup(OCPval.G,OCPval.F,A,l,u,'eps_abs', 1e-04, 'eps_rel', 1e-04,'polish','true');
end

