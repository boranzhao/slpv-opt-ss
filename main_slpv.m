%% Desing of switching LPV (SLPV) controller (with improved local performance)

% Idea of the method:
% (1) For LPV plants with affine parameter-dependence, multi-convexity (MC)
% concept is used to get a finite number of LMIs using affine Lyapunov variables.
% (2) Different cost functions can be selected in SLPV controller design: 
%     1) max{gamma^(j)} for traditional SLPV controllers design
%     2) weighted sum of gamma^(j) for SLPV controllers with improved local performance


Initialization
if design_para.XY_PD ~=3
    design_para.ss_improve =0;
    design_para.max_iter = 1; % no need to do iterations on ss condition
end

switch plant_index
    case 5
        tau1_vec = -0.9:0.05:0.9; step = 0.1;
    case {1 2} 
        tau1_vec = -0.7:0.1:0.7; step = 0.1;
    case 8
        tau1_vec = 400:50:1000; step = 50; 
    case 6
        tau1_vec = -2.8:0.4:2.8; step = 0.4; 
    case 10 
        tau1_vec = -1.6:0.2:1.6; step = 0.2; 
end
if design_para.full_search == 1    
    if plant_index ==12
        [tau1_1,tau1_2] = meshgrid(1.464:0.02:13.642, 1.464:0.02:13.642);
    else
        [tau1_1,tau1_2] = meshgrid(-0.99:0.01:0.99, -0.99:0.01:0.99);
    end
%     tau1_1 = tau1_1(tau1_1<tau1_2);
%     tau2_2 = tau1_1(tau1_1<tau1_2);
    GamM = zeros(size(tau1_1)); 
    gam_full = ones(numel(tau1_1),regnum+4)*NaN;
    id = 0;
    tic;
    for i = 1: size(tau1_1,1)
        for j= 1:size(tau1_1,2)
            if tau1_1(i,j) >= tau1_2(i,j)
                GamM(i,j) = NaN;
            else      
                id = id + 1;
                tau1 = [tau1_1(i,j) tau1_2(i,j)];
                theta_info.Theta1 = [Theta1(1,1) tau1_2(i,j);tau1_1(i,j) Theta1(2,2)];
                if design_para.LPVdesign_method == 1
                    [OptRst,sol] = SLPV_Basic_Affine_MC(Gasym,Ga_ctmat,theta_info,design_para,opt_lmilab);
                elseif  design_para.LPVdesign_method == 2                    
                    [OptRst,sol] = SLPV_Proj_Affine_MC(Gasym,Ga_ctmat,theta_info,design_para,opt_lmilab);
                end
                % note that all the variables with postfix _k are the value of associated variables after iteration k 
%                 disp(['Gam: ' num2str(OptRst.Gam) '    gam:' num2str(OptRst.gam)  '    tau1:' num2str(tau1)]);
                GamM(i,j) = OptRst.Gam;
                gam_full(id,1) = OptRst.Gam;
                gam_full(id,2:2+regnum-1) = OptRst.gam;
                gam_full(id,2+regnum) = OptRst.obj;
                gam_full(id,2+regnum+1:end)= tau1;
            end
        end
    end 
    time_full = toc; % full-search time
    gam_full = gam_full(~isnan(gam_full(:,1)),:);
    [Gam_min,id_min] = min(gam_full(:,1));
    disp('Minimal \gamma was obtained as');
    disp(['Gam: ' num2str(gam_full(id_min,1)) '    gam:' num2str(gam_full(id_min,2:2+regnum-1))  '    tau1:' num2str(gam_full(id_min,2+regnum+1:end))]);
    return;      
%     tau1_full = ones(1000,2)*inf;
%     i = 1;
%     for tau1_1 =  tau1_vec
%         for tau1_2 = tau1_1+step:step:tau1_vec(end) % (tau1_1/10+0.1:0.1:0.9)*10
%             tau1_full(i,:) = [tau1_1 tau1_2];
%             i = i + 1;
%         end
%     end
%     id = find(tau1_full(:,1)~=inf);
% %     tau1_full = tau1_full(1:sum(tau1_full(:,1)~=0 | tau1_full(:,2)~= 0),:);
% %     tau1_full = [1.9937 2.6424];
%     tau1_full = tau1_full(id,:);
%     gam_full = zeros(size(tau1_full,1),regnum+4);
%     % gam obj tau1 
%     for i = 1: size(tau1_full,1)
%         tau1 =  tau1_full(i,:);
%         theta_info.Theta1 = [Theta1(1,1) tau1(2);tau1(1) Theta1(2,2)];
%         [OptRst,sol] = Copy_of_SLPV_Proj_Affine_MC(Gasym,Ga_ctmat,theta_info,design_para,opt_lmilab);
%         % note that all the variables with postfix _k are the value of associated variables after iteration k 
%         disp(['Gam: ' num2str(OptRst.Gam) '    gam:' num2str(OptRst.gam)  '    tau1:' num2str(tau1)]);
%         gam_full(i,1) = OptRst.Gam;
%         gam_full(i,2:2+regnum-1) = OptRst.gam;
%         gam_full(i,2+regnum) = OptRst.obj;
%         gam_full(i,2+regnum+1:end)= tau1;
%     end 
%     time_full = toc; % full-search time
%     return;    
else
    if design_para.ss_improve == 0
        design_para.max_iter = 1;
    end
    X_k =[]; Y_k = []; Gam_k =[];
    OptRst_k = [];
    gam_iter = ones(design_para.max_iter,regnum+1)*NaN;
    tic;
    for k=1:design_para.max_iter
        if design_para.LPVdesign_method == 1
            [OptRst,sol] = SLPV_Basic_Affine_MC(Gasym,Ga_ctmat,theta_info,design_para,opt_lmilab,OptRst_k);
        elseif  design_para.LPVdesign_method == 2                    
            [OptRst,sol] = SLPV_Proj_Affine_MC(Gasym,Ga_ctmat,theta_info,design_para,opt_lmilab,OptRst_k);
        end
        disp(['Iter: ' num2str(k) '   Gam: ' num2str(OptRst.Gam) '    gam:' num2str(OptRst.gam)])
        gam_iter(k,1) = OptRst.Gam;
        gam_iter(k,2:end) = OptRst.gam;
        imprvchk = k>5 && abs(gam_iter(k-5,1)-gam_iter(k,1))/gam_iter(k-5,1) < design_para.tol;
        if imprvchk
            disp('SS iteration terminated due to no significant improvement in Gamma vlaue!');
            break;
        end
        X_k = OptRst.X;Y_k = OptRst.Y;  
        OptRst_k = OptRst; 
    end  
    gam_iter = gam_iter(~isnan(gam_iter(:,1)),:);
    toc;
end
if design_para.SwLogic == 0 % && gam_min_method == 2
    Gam_NS = OptRst.Gam;
end  
if design_para.gam_min_method == 1
    gam_iter_max = gam_iter;
elseif design_para.gam_min_method == 2
    gam_iter_ws = gam_iter;
end
% return;
%% check whether original conditons for controller design are satisfied 
% to use the result from LPV controller for validation on switching LPV
% controller 
% OptRst1 = OptRst;
% design_para.SwLogic = 1;
% theta_info.Theta1 = [1.464 8;7.5 13.642]; 
% design_para.REGID = [1 1 1;2 2 1];
% Xtmp(:,:,:,1) = OptRst.X;
% Xtmp(:,:,:,2) = OptRst.X;
% OptRst1.X = Xtmp;
% Ytmp(:,:,1,1) = OptRst.Y;
% Ytmp(:,:,1,2) = OptRst.Y;
% OptRst1.Y = Ytmp;
% OptRst1.lam = [OptRst.lam OptRst.lam];
% OptRst1.sig = [OptRst.sig OptRst.sig];
% OptRst1.gam = [OptRst.gam OptRst.gam];
% OptRst1.mu = [OptRst.mu OptRst.mu];
%%
% [consts_eval,satisfy,ineq_num] = SLPV_ss_check_MC(Gasym,Ga_ctmat,theta_info, design_para, OptRst);
% satisfy      

%% Run controller Construction
ControllerConstruction

Gam = max(OptRst.gam);
% filter parameter
swtch = 1;
alpha = 0.1;
