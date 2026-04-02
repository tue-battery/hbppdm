function [wp,lambda,flag,k,rho_e,rho_p] = HBPPDM_var_reset(OCPval,epsilon,niter,alpha,gamma,wp,lambda)  
    %initialize commonly used inverses and multiplications
    invGAT = OCPval.iG*(OCPval.AA.');
    invGF = OCPval.iG*OCPval.F;
    AinvGAT = OCPval.AA*invGAT;
    
    %initialize w and lambda with zeros.
    if isempty(lambda)
        lambda  = -alpha*(AinvGAT\(OCPval.bb-OCPval.AA*invGF));
    end
    
    if isempty(wp)
        wp = zeros(length(OCPval.AA),1);
    end

    lambdam = zeros(length(OCPval.bb),1);
    

    n = 1;
    
    for k = 1:niter
       b = n./(n+gamma);
       wp = max(OCPval.wlb,min(OCPval.wub,-(invGAT*lambda+invGF)));  
       z_k = OCPval.AA*wp-OCPval.bb;
       
       %check convergence
       rho_p(k) = norm(z_k);
       if rho_p(k) <= epsilon
           break
       end
       
       % Eq. 16c \w phi22 = A*inv(G)*A'
%        a(k) = (k-1)/(k+alpha);
       lambdah = lambda +b*(lambda-lambdam);
       lambdam = lambda;
       lambda = lambdah + alpha*((AinvGAT)\(z_k));        
       
       if (lambdah-lambda)'*(lambda-lambdam) >= 0
           n = -1;
       end
       
       % Eq. 16
          
       n = n+1;
    end
    
    rho_e = norm(z_k);

    %check if solution is found and set output
    if k == niter
        flag = 1;
    else
        flag = 0;
    end

end