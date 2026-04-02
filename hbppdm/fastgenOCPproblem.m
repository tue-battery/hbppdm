function [OCPval] = fastgenOCPproblem(A,B,Q,P,R,N)
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
        
    iQ = inv(Q);
    iP = inv(P);
    iR = inv(R);

    %% G & iG
%     OCPval.G = sparse((N+1)*n+N*m);
%     OCPval.iG = sparse((N+1)*n+N*m);
%     for var_for = 1:N
%        OCPval.G(1+(var_for-1)*n:var_for*n,1+(var_for-1)*n:var_for*n) = Q;
%        OCPval.iG(1+(var_for-1)*n:var_for*n,1+(var_for-1)*n:var_for*n) = iQ;
%     end
%     OCPval.G(n*N+1:n*(N+1),n*N+1:n*(N+1)) = P;
%     OCPval.iG(n*N+1:n*(N+1),n*N+1:n*(N+1)) = iP;
%     for var_for = 1:N
%        OCPval.G((N+1)*n+1+(var_for-1)*m:(N+1)*n+var_for*m,(N+1)*n+1+(var_for-1)*m:(N+1)*n+var_for*m) = R;
%        OCPval.iG((N+1)*n+1+(var_for-1)*m:(N+1)*n+var_for*m,(N+1)*n+1+(var_for-1)*m:(N+1)*n+var_for*m) = iR;
%     end
    
    q = zeros((N+1)*n+N*m,1);
    iq = zeros((N+1)*n+N*m,1);
    Q_diag = diag(Q);
    iQ_diag= diag(iQ);
    R_diag = diag(R);
    iR_diag= diag(iR);
    for var_for = 1:N
        q(1+(var_for-1)*n:var_for*n,1) = Q_diag;
        iq(1+(var_for-1)*n:var_for*n,1) = iQ_diag;
    end
    q(n*N+1:n*(N+1),1) = diag(P);
    iq(n*N+1:n*(N+1),1) = diag(iP);
     for var_for = 1:N
        q((N+1)*n+1+(var_for-1)*m:(N+1)*n+var_for*m,1) = R_diag;
        iq((N+1)*n+1+(var_for-1)*m:(N+1)*n+var_for*m,1) = iR_diag;
     end

    OCPval.G = sparse(1:(N+1)*n+N*m,1:(N+1)*n+N*m,q);
    OCPval.iG = sparse(1:(N+1)*n+N*m,1:(N+1)*n+N*m,iq);
    clear q iq;
    %%   
%     gamma_x = -speye((N+1)*n);
%     gamma_x(1:n,1:n) = eye(n);
%     
%     for k = 1:N
%         gamma_x(k*n+1:k*n+n,(k-1)*n+1:k*n) = A;
%     end
    
    [A_row,A_col,A_var] = find(A);
    [B_row,B_col,B_var] = find(B);
    len_A = length(A_var);
    len_B = length(B_var);
    for k = 1:N
        gamma_row(1+(k-1)*len_A:k*len_A,1) = A_row +k*n;
        gamma_col(1+(k-1)*len_A:k*len_A,1) = A_col +(k-1)*n;
        gamma_val(1+(k-1)*len_A:k*len_A,1) = A_var;
        gamma_row(1+N*len_A+(k-1)*len_B:N*len_A+k*len_B,1)=B_row + k*n;
        gamma_col(1+N*len_A+(k-1)*len_B:N*len_A+k*len_B,1)=B_col + (k-1)*m + n*(N+1);
        gamma_val(1+N*len_A+(k-1)*len_B:N*len_A+k*len_B,1)=B_var;
    end
    gamma_row = [gamma_row; (1:n*(N+1))'];
    gamma_col = [gamma_col; (1:n*(N+1))'];
    gamma_val = [gamma_val; ones(n,1); -ones(n*N,1)];
    
    OCPval.gamma_row = gamma_row;
    OCPval.gamma_col = gamma_col;
    OCPval.gamma_val = gamma_val;
    OCPval.AA = sparse(gamma_row,gamma_col,gamma_val,OCPval.nStates*(N+1),OCPval.nStates*(N+1)+OCPval.nInputs*N);

end