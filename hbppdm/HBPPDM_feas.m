function [wp,lambda,flag,k] = HBPPDM_feas(OCPval,epsilon,niter,alpha,gamma)
    %initialize commonly used inverses and multiplications
    invGAT = OCPval.iG*(OCPval.AA.');
    invGF = OCPval.iG*OCPval.F;
    AinvGAT = OCPval.AA*invGAT;
    %initialize w and lambda with zeros.
    lambda  = -AinvGAT\(OCPval.bb-OCPval.AA*invGF);
    lambdam = zeros(length(OCPval.bb),1);
    lambda_zeros = zeros(length(OCPval.AA),1);
    wp = zeros(length(OCPval.AA),1);
    n = 1;
    for k = 1:niter
      if OCPval.bb'*lambda + OCPval.wub'*max(lambda_zeros,-OCPval.AA'*lambda) + OCPval.wlb'*min(lambda_zeros,-OCPval.AA'*lambda) < 0
           printf('unfeasible')
           break;
       end
       b = n./(n+gamma);
       wp = max(OCPval.wlb,min(OCPval.wub,-(invGAT*lambda+invGF)));  
       z_k = OCPval.AA*wp-OCPval.bb;
       
       %check convergence
       
       if norm(z_k) <= epsilon
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
    %check if solution is found and set output
    if k == niter
        flag = 1;
    else
        flag = 0;
    end
end