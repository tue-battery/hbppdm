%% HBPPDM EXAMPLES
% three cases:
% - pendulum input constraint           (AUTOMATICA)
% - pendulum input and state constraint (AUTOMATICA)
% - AFTI Reference tracking problem     (AUTOMATICA)
% - eco-driving                         (IFAC WC 2026)
% - eco-driving (MAASTRICHT ROUTE)      (IFAC WC 2026)

clear;
close all;

addpath("ecodriving_maastricht\");
addpath("hbppdm\")
addpath("ocp_problems\")

%% SETTINGS
%Case:
% 1 - pendulum i constraint, 
% 2 - pendulum i&s constraint, 
% 3 - ATFI-16
% 4 - eco-driving
% 5 - eco-driving Maastricht route
benchmark_no = 5;

%Horizon length
N = 1000;

%Maximum iterations hbppdm
max_iter = 1e5;

%Final residual for hbppdm
resid = 1e-4;

% also run osqp
run_osqp = 0;

% also run quadprog
run_qp = 0;

% also run cplex
run_cplex = 0;

% also run scs
run_scs = 0;

%% CASES 1 - 2 - 3
if benchmark_no == 1 || benchmark_no == 2 || benchmark_no == 3

%select case
switch benchmark_no
    case 1
        OCPval = pendulumunc(N);
        bcase = "pendulumi";
        fprintf('GENERATING PENDULUM (input constraint)\n')
        alpha = 1.96;
        gamma = 28.22;
    case 2
        OCPval = pendulumc(N);
        bcase = "pendulumis";
        fprintf('GENERATING PENDULUM (input & state constraint)\n')
        alpha = 1.96;
        gamma = 28.22;
    case 3
        OCPval = ATFI16func(N);
        bcase = "ATFI16";
        fprintf('GENERATING ATFI (input & state constraint & reference)\n')
        alpha = 1.96;
        gamma = 28.22;
end
    
    fprintf('RUNNING OCP WITH HORIZON: %2.0f. CASE: ' + bcase + ' \n', N);

    %% HBPPDM    
    tic
    [w_hbppdm,lambda_hbppdm,flag,iter,~,~] = HBPPDM_var_reset(OCPval,resid,max_iter,alpha,gamma,[],[]);
    time_hbppdm = toc;
    pf.hbppdm = norm(OCPval.AA*w_hbppdm-OCPval.bb);
    cost_hbppdm = 0.5*w_hbppdm'*OCPval.G*w_hbppdm + w_hbppdm'*OCPval.F;
    
    %% Quadprog
    if run_qp
    tic
    [w_qp,fval,exitflag,output] = quadprog(OCPval.G,OCPval.F,[],[],OCPval.AA,OCPval.bb,OCPval.wlb,OCPval.wub);
    time_qpt = toc;
    pf.qp = norm(OCPval.AA*w_qp-OCPval.bb);
    cost_qp = 0.5*w_qp'*OCPval.G*w_qp + w_qp'*OCPval.F;
    end
    
    %% cplex
    if run_cplex
    tic;
    w_cpl = cplexqp(OCPval.G,OCPval.F,[],[],OCPval.AA,OCPval.bb,OCPval.wlb,OCPval.wub);
    time_cpl = toc;
    pf.cpl = norm(OCPval.AA*w_cpl-OCPval.bb);
    cost_cpl = 0.5*w_cpl'*OCPval.G*w_cpl + w_cpl'*OCPval.F;
    end

    %% OSQP
    if run_osqp
    m = prepareOSQP(OCPval);
    tic
    results = m.solve();
    time_osqp = toc;
    w_osqp = results.x;
    pf.osqp = norm(OCPval.AA*w_osqp-OCPval.bb);
    cost_osqp = 0.5*w_osqp'*OCPval.G*w_osqp + w_osqp'*OCPval.F;
    end
        
    %% SCS
    if run_scs
    [scsdata,scscone] = prepareSCS(OCPval);
    ScsSettings = struct;
    
    tic
    [x, y, s, info] = scs(scsdata, scscone, []);
    w_scs = x;
    time_scs = toc;
    pf.scs = norm(OCPval.AA*w_scs-OCPval.bb);
    cost_scs = 0.5*w_scs'*OCPval.G*w_scs + w_scs'*OCPval.F;
    end
    
end


%% ECO-DRIVING CASES

if benchmark_no == 4
    iter = 50;
    wp2 = [1; repmat([2.5;.1; sqrt((0.5*0.15)/2.5) ;0],N,1)];
    OCPval = [];
    lambda = [];
    
    % SQP layer
    for k_hbppdm = 1:iter
    
    %linearize problem
    [OCPval] = ecodrivingOCP(OCPval, wp2, N);
    OCPval.wp = wp2;
    
    % run hbppdm
    tic
    [wpp,lambda,flag,hbppdm_iter(k_hbppdm),rho_sqp, rho_p{k_hbppdm}] = HBPPDM_var_reset(OCPval,1e-4,100,1.95,15,wp2,lambda);
    t_sqp_hbppdm(k_hbppdm) = toc;
            

    epsilon = norm(wpp-wp2);
    wp2 = wpp;
    
    % save for figures
    rho_hbppdm(k_hbppdm) = norm(OCPval.AA*wpp-OCPval.bb);
    epsilon_hbppdm(k_hbppdm) = epsilon; 
    
    if epsilon < 1e-4
        break
    end
    end

end %end case 4


if benchmark_no == 5
    iter = 50;
    Ts = 20;
    [mstr_data] = maastricht_load(Ts);
    N = mstr_data.N;

    wp2 = [.13; repmat([2.0;.1; sqrt((0.5*0.13)/2.5) ;0],N,1)];
    OCPval = [];
    lambda = [];
    
    % SQP layer
    for k_hbppdm = 1:iter
    
    %linearize problem
    [OCPval] = ecodrivingOCP_maastricht(OCPval, wp2, mstr_data);
    OCPval.wp = wp2;
    
    % run hbppdm
    tic
    [wpp,lambda,flag,hbppdm_iter(k_hbppdm),rho_sqp, rho_p{k_hbppdm}] = HBPPDM_var_reset(OCPval,1e-4,100,1.95,15,wp2,lambda);
    t_sqp_hbppdm(k_hbppdm) = toc;

    
    epsilon = norm(wpp-wp2);
    wp2 = wpp;
    
    % save for figures
    rho_hbppdm(k_hbppdm) = norm(OCPval.AA*wpp-OCPval.bb);
    epsilon_hbppdm(k_hbppdm) = epsilon; 
    
    if epsilon < 1e-4
        break
    end
    end

end %end case 5
