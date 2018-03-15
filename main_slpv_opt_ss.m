%% Desing of switching LPV (SLPV) controller with (improved local performance and) optimized switching surfaces

% Idea of the method:
% (1) For LPV plants with affine parameter-dependence, multi-convexity (MC)
% concept is used to get a finite number of LMIs using affine Lyapunov variables.
% (2) Different cost functions can be selected in SLPV controller design: 
%     1) max{gamma^(j)} for traditional SLPV controllers design
%     2) weighted sum of gamma^(j) for SLPV controllers with improved local performance%  (3)  An algorithm based on particle swarm optimization is used to optimize the SSs.                     

% clc;
% clear;
Initialization
ga_pso_sel = 2;     % 1 for using genetic algorithm; 2 for using particle swarm optimization (PSO)
%% constraints on SSVs
nvars = (regnum1-1)*2+(regnum2-1)*2;
tau1_len = (regnum1-1)*2; tau2_len = nvars-tau1_len;
if tau2_len ==0
    total_const_num = tau1_len-1 +2;
else
    total_const_num = tau1_len-1 +2 + tau2_len-1 +2;
end
Aineq = zeros(100,nvars);
bineq = zeros(100,1);
for index = 1:tau1_len-1
    Aineq(index, [index index+1]) = [1 -1];
    bineq(index,1) = -ss_dist_min(1);
end

for index = 1:max(0,tau2_len-1)
    Aineq(tau1_len-1+index,[tau1_len+index tau1_len+1+index]) = [1 -1];
    bineq(tau1_len-1+index,1) = -ss_dist_min(2);
end
cur_const_num = tau1_len-1+max(0,tau2_len-1) +1;
Aineq(cur_const_num,1) = -1; bineq(cur_const_num,1) = -(theta_info.Theta1(1,1)+ss_dist_min(1));
cur_const_num = cur_const_num +1;
Aineq(cur_const_num,tau1_len) = 1; bineq(cur_const_num,1) =  theta_info.Theta1(end,end)-ss_dist_min(1);

if tau2_len > 0 
    cur_const_num = cur_const_num +1;
    Aineq(cur_const_num,tau1_len+1) = -1; bineq(cur_const_num,1) = -(theta_info.Theta2(1,1)+ss_dist_min(2));
    cur_const_num = cur_const_num +1;
    Aineq(cur_const_num,nvars) = 1; bineq(cur_const_num,1) =  theta_info.Theta2(end,end)-ss_dist_min(2);    
end

% %% for engine example, limit the constraint of tau_1,2, tau_1,3
% cur_const_num = cur_const_num+1;
% Aineq(cur_const_num,2)= 1; bineq(cur_const_num) = 2.5;
% cur_const_num = cur_const_num+1;
% Aineq(cur_const_num,3)= -1; bineq(cur_const_num) = -10;

% lower and upper bounds
LB = [theta_info.Theta1(1,1)*ones((regnum1-1)*2,1); theta_info. Theta2(1,1)*ones((regnum2-1)*2,1)];
UB = [theta_info.Theta1(end,end)*ones((regnum1-1)*2,1); theta_info. Theta2(end,end)*ones((regnum2-1)*2,1)];

% pso parameters 
kappa = 1;
phi1 = 2.05;
phi2 = 2.05;
phi = phi1+phi2;
chi = 2*kappa/abs(2-phi-sqrt(phi^2-4*phi));
c1 = chi*phi1;             % Cognitive Attraction: same as matlab default
c2 = chi*phi2;             % Social Attraction: same as matlab default              
nPop = 10;%nvars*10;           % Population size (Swarm Size), matlab default:  min(100,10*nvars)
c0 =  chi;                 % Intertia coefficient
MaxIter = 50;               % Generations

tic;
cost_fun  = @(x) gamma_cal_slpv(x, Gasym, Ga_ctmat, theta_info, design_para,opt_lmilab);
if ga_pso_sel == 1
    [tau_opt,Gam_opt,exitflag,output] = ga(cost_fun,nvars,Aineq,bineq,[],[],LB,UB,[]);%
elseif ga_pso_sel == 2
    options = psooptimset('TolFun',5e-4,'TolCon',1e-5,'CognitiveAttraction',c1,'SocialAttraction',c2,...
        'Generations',MaxIter,'PopulationSize',nPop,'StallGenLimit',8,'PlotFcns',{@psoplotbestf,@psoplotswarm}); 
    [tau_opt,Gam_opt,exitflag,output,population,scores,iterations] = pso_ssopt(cost_fun,nvars,Aineq,bineq,[],[],LB,UB,[],options);%
end
time_ga_pso = toc;
%% plot the iteration result
close all;
figure;
subplot(2,1,1)
plot(iterations.fGlobalBest,'r-o','linewidth',1.5); hold on;
plot(max(iterations.fLocalBest),'k-.*','linewidth',1.5);
plot(mean(iterations.fLocalBest),'b--*','linewidth',1.5);
xlabel('Iteration');
ylabel('Cost function');
legend('global best','worst of local bests','mean of local bests');
% goodplot
subplot(2,1,2)
plot(squeeze(iterations.xGlobalBest)','Linewidth',1.5); 
xlabel('Iteration');
ylabel('SSVs ');
legend({'\tau_{1,1}','\tau_{1,2}','\tau_{1,3}','\tau_{1,4}'},'orientation','hori');

%%
% goodplot([7 6])
% print -painters -dpdf -r150 pso_iter.pdf


