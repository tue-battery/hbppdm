function [wp,lambda,cost,iter] = HBPPDM_blksqp(H,iH,OCPval,opts)
%Wrapper for optimal control problems with blockdiagonal costfunction.
%Yannick Heuts, 12-april-2021

% Unpack OCPval.
G = OCPval.G;
F = OCPval.F;
AA= OCPval.AA;
bb= OCPval.bb;
wlb=OCPval.wlb;
wub=OCPval.wub;

% Setting options
nsqp = opts.nsqp;
nhbppdm = opts.nhbppdm;
alpha = opts.alpha;
b = opts.b;
eps1 = opts.eps1;
eps2 = opts.eps2;

% initializing global values
Phi1 = iH*(AA.');
Phi2 = AA*Phi1;
lambdam = zeros(length(bb),1);
%wp = [OCPval.xinit;zeros(OCPval.N*(OCPval.nStates+OCPval.nInputs),1)];
wp = OCPval.wp;
cost = 0;
pf = inf;
% sqp iteration
for j = 1:nsqp
    F2 = (G-H)*wp+F;
    Phi3 = iH*F2;
    if j == 1
        lambda  = -Phi2\(bb-AA*Phi3);
    end
    % Heavy ball projected primal dual method
    [wp,lambda,iter(j,1),pf] = HBPPDM_(Phi1,Phi2,Phi3,AA,bb,wlb,wub,alpha,b,nhbppdm,min(pf,eps1),lambda);
    % Checking cost
    cost(j+1) = 0.5*wp'*G*wp+wp'*OCPval.F;
    % If cost is within tolerance, solve one more subproblem with higher
    % tolerance for primal feasibility
    if abs((cost(j) - cost(j+1))/cost(j)) < 1e-4
        F2 = (G-H)*wp+F;
        Phi3 = iH*F2;
        lambda  = -Phi2\(bb-AA*Phi3);
        [wp,lambda,iter(j,1)] = HBPPDM_(Phi1,Phi2,Phi3,AA,bb,wlb,wub,alpha,b,nhbppdm,min(pf,eps2),lambda);
        break
    end   

end
end

function [wp,lambda,k,pf] = HBPPDM_(Phi1,Phi2,Phi3,AA,bb,wlb,wub,alpha,b,niter,epsilon,lambda)
    %Heavy Ball Projected Primal Dual Method
    %Yannick Heuts, 12-april-2021
    n = 1;
    lambdam = zeros(length(lambda),1);
    for k = 1:niter
       beta = n./(n+b);
       wp = max(wlb,min(wub,-(Phi1*lambda+Phi3)));  
       z_k = AA*wp-bb;
       %check convergence
       pf = norm(z_k);
       if pf <= epsilon
           break
       end
       if mod(k,100) == 1
           disp("rho: " + string(pf))
       end
       lambdah = lambda +beta*(lambda-lambdam);
       lambdam = lambda;
       lambda = lambdah + alpha*(Phi2)\(z_k);               
       if (lambdah-lambda)'*(lambda-lambdam) >= 0
           n = -1;
       end
       n = n+1;
    end
end

