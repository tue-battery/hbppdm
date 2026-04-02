function [OCPval] = genOCPproblem(ssd,Q,P,R,xub,xlb,uub,ulb,N,xinit)
% ssd - discrete state space object of plant
% Q - penalty matrix for state
% P - penalty matrix for final state
% R - penalty matrix for input
% xub - upper bound for state (N*n x1 double)
% xlb - lower bound for state (N*n x1 double)
% uub - upper bound for input (1x1 double)
% uub - lower bound for input (1x1 double)
% N - horizon
% xinit - initial condition 



% [OCPval.G, OCPval.AA, OCPval.bb] = gen_OCP_matrices(ssd.A,ssd.B,Q,R,P,N,xinit);
%%
    n = length(Q);
    m = length(R);
    
    OCPval.nStates =length(Q);
    OCPval.nInputs =length(R);
    OCPval.N = N;
    
    Plen = length(P);
    
    iQ = inv(Q);
    iP = inv(P);
    iR = inv(R);
    % Building G
    OCPval.G = zeros(N*length(Q)+length(P)+N*length(R));
    OCPval.G(1:n*N,1:n*N) = kron(eye(N),Q);
    OCPval.G(n*N+1:n*N+Plen,n*N+1:n*N+Plen) = P;
    OCPval.G(n*N+Plen+1:n*N+Plen+m*N,n*N+Plen+1:n*N+Plen+m*N) = kron(eye(N),R);
    
    OCPval.iG = zeros(N*length(Q)+length(P)+N*length(R));
    OCPval.iG(1:n*N,1:n*N) = kron(eye(N),iQ);
    OCPval.iG(n*N+1:n*N+Plen,n*N+1:n*N+Plen) = iP;
    OCPval.iG(n*N+Plen+1:n*N+Plen+m*N,n*N+Plen+1:n*N+Plen+m*N) = kron(eye(N),iR);
    
    gamma_x = eye((N+1)*n);
    gamma_x(n+1:end,n+1:end) = -gamma_x(n+1:end,n+1:end);
    
    for k = 1:N
        gamma_x(k*n+1:k*n+n,(k-1)*n+1:k*n) = ssd.A;
    end
    
    gamma_u = [zeros(n,m*(N));kron(eye(N),ssd.B)];
    OCPval.AA = [gamma_x, gamma_u];
    OCPval.bb = [xinit; zeros(n*N,1)];
    
    
%%

    for k = 0:length(xub)-1
       OCPval.xub(1+n*k:n+n*k,1) = xub(k+1,:)';
       OCPval.xlb(1+n*k:n+n*k,1) = xlb(k+1,:)';
    end

    OCPval.uub = kron(ones(N,1),uub);
    OCPval.ulb = kron(ones(N,1),ulb);
    OCPval.wub = [OCPval.xub; OCPval.uub];
    OCPval.wlb = [OCPval.xlb; OCPval.ulb];

    OCPval.F = zeros(length(OCPval.G),1);

    
%%
    OCPval.G = sparse(OCPval.G);
    OCPval.iG= sparse(OCPval.iG);
    OCPval.AA= sparse(OCPval.AA);
    
    OCPval.xinit = xinit;
    
    OCPval.A = ssd.A;
    OCPval.B = ssd.B;
    OCPval.C = ssd.C;

end