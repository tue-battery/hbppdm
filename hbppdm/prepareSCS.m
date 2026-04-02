function [data,Cone] = prepareSCS(OCPval)
% build required matrices
data.P = OCPval.G;
data.A = [OCPval.AA; sparse(1,(OCPval.N+1)*OCPval.nStates+OCPval.N*OCPval.nInputs); -speye((OCPval.N+1)*OCPval.nStates+OCPval.N*OCPval.nInputs)];
% data.A = [OCPval.AA];
data.c = OCPval.F;
data.b = [OCPval.bb;1;zeros((OCPval.N+1)*OCPval.nStates+OCPval.N*OCPval.nInputs,1)];
% data.b = [OCPval.bb];

Cone.l = 0;
Cone.z = (OCPval.N+1)*OCPval.nStates;
Cone.bu = OCPval.wub(1:end);
Cone.bl = OCPval.wlb(1:end);
Cone.bsize = length(OCPval.wub);
end



