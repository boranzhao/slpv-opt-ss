%% parameter initialization                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               global opt_lmilab Ga_ctmat design_para 
% public parameters
design_para.solver = 1;             % 1 for Matlab LMI lab, 2 for yalmip (currently not available)
design_para.XY_PD = 3;              % 0 for parameter-independent X and Y, 1 for parameter-dependent X and constant Y; 2 for parameter-dependent Y and constant X, 3 for parameter-dependent X and Y
design_para.full_search = 0;        % 0 for nothing, 1 for full search, 2 for GA/PSO

design_para.LPVdesign_method = 2;   %1 for basic characterization, 2 for projected characterization. See [Apkarian 1998] for details
design_para.gam_min_method = 1;     % 1 for minimization of max gam_j ; 2 for minimization of weighted sum of gam_j based on subset sizes
design_para.SwLogic = 1;            % 1 for non-switching, 1 for hesteresis switching
design_para.ss_opt_part = 1;        % just optimize part of switching surfaces 

design_para.ss_improve = 0;         % flag for determining whether to use convex-concave decomposition to improve ss condition. ?not included in the paper)
design_para.ss_cond = 2;            % only used for XY_PD=3, 1 for Y=Y1, X>=X1; 2 for X=X1 Y<=Y1;
sign_para.max_iter = 1;             % maximum iteration in using convex-concave decomposition to improve ss condition. 

design_para.num_cost_gain = 5e-7*[1 1 1]; % for penalizing the magnitude of X&Y M&N and sig
design_para.overlap_gain = 1e-3/5200*0;   % for penalizing overlapped area between adjacent subsets

design_para.sig_limit = 1e6;        % for avoiding too large sig that may cause numerical issues
design_para.XY_norm_limit = 1e5;    % for avoiding too large X or Y that may cause numerical issues
design_para.lam_mu_limit = 1e6;     % for avoiding too large mu that may cause numerical issues
design_para.overlap_limit = 5200*5; % limit of overlapped region between adjacent subsets

design_para.eps_div = 1e-9;         % smaller number to avoid division by zero. 
design_para.eps_eq = 1e-6;          % for imposing equality between two matrices using LMIs, 2-norm of the difference between the two matrices will be smaller than this value 
design_para.LMI_relax = [design_para.eps_eq design_para.eps_eq]*2; %LMI_relax(1) used for initial value computation, LMI_relax(2) used for iteration 
design_para.tol = 1e-5;             % relative tolerance for terminating iteration to improve SS condtion

% multiconvexity
design_para.mc_relax = 1;           % 0: gridding technique; 1: relax the muliti-convexity constraints by using  -(lambda_0 +
%  sum_theta_i^2*lambda_i+sum_theta_i*lambda_i); 2: no relax: i.e. using 0

% method for convexifying originally nonconvex switching surface conditions, 
% 1 for using same lmivar to impose equaility of either X or Y at the SS;
% 2 for using inequality constraints to impose equaility of either X or Y at the SS;
% 3 for adding an extra constraint inv(Y)>=eps*I
design_para.ss_cvx.method = 1;
design_para.ss_cvx.eps = 1e-7; 

%% settings for the LMI Lab solver
opt_lmilab(1)= 1e-5;    % relative accuary on the optimal value
opt_lmilab(2)= 600;     % Number of Iteration
opt_lmilab(4) = 10;     % J, the code terminates when the objective has not decreased by more than the desired relative accuracy during the last J iterations
opt_lmilab(5)=  1;   

%% Plant selection 
% 1 for numerical example, 
% 2 for engine control example
%%  numerical example
plant_index = 1;
Theta1 = [-1 0.1;-0.1 1];
% Theta1 = [-1 -0.3;-0.4 0.3;0.2 1];
% % Theta1 = [-1 0.810;0.800 1];          % optimal SS 
% Theta1 = [-1 0.8203;0.8101 1];          % optimal SS for gamma_bar equal 1.02*traditional SLPV
% Theta1 = [-1 0.81214;0.80096 1];        % optimal SS even for traditonal SLPV controller

%% engine control example
% 1 subset
% Theta1 = [1.464  13.642];

% 2 subsets 
% Theta1 = [1.464 5.2947;2.4202 13.642]; % pso optimization

% 3 subsets
% Theta1 = [1.464 6;5.5 10;9.5 13.642]; % heuristically: used in CCTA paper
% Theta1 = [1.464 6.5;6 10.5;9.5 13.642]; % heuri


% 4 subsets
%  Theta1 = [1.464 5;4 8;7  11;10 13.642]; 

% pso
% Theta1 = [1.464  4.599;3.1219 10.469;4.8080 13.642]; % pso,XY:3, proj
% Theta1 = [1.464  4.7702;3.04243  10.7685;5.01555  13.642]; % pso,XY:3, proj
% Theta1 = [1.464  4.631;3.121 10.397;4.795 13.642]; % pso,XY:3, proj, pso_3subs_XY3_Proj_2perc3
%  Theta1 = [1.464 6.5;4.55057 10.5554;8.45505 13.642];  %  pso,XY:3, proj, pso_3subs_XY3_Proj_1perc
%  Theta1 = [1.464 1.75883 ;1.64252 11.553;11.1491 13.642]; 
%  Theta1  = [1.464 6.651;4.551  10.556;8.455 13.642]; % pso_3subs_XY3_max
% Theta1  = [1.464 4.3374;1.8779  9.9162;4.7507 13.642]; %pso_3subs_XY3_Proj_max1

% for test
% Theta1 = [6.1 6.2];
gs_num = 1; 

[Gasym,theta_fcn_info,gs_num,theta_info,num_gain] = plant_sel(plant_index,Theta1,gs_num); % use Plant_sel.m for automotive engine example 
ss_dist_min = theta_info.ss_dist_min;% the minimum distance between any two adjacent switching surfaces 
theta_info.agg_subset = ones(1,size(theta_info.Theta1,1)*size(theta_info.Theta2,1)); % denote which subset is considered as aggressive subsets
design_para.subset_weight = [];%[1 0 0]; % will use this weight if not empty [0 1 0]; %[0.2 0.6 0.2];
if ~exist('Gam_NS') || design_para.SwLogic ==0
    Gam_NS = [];
    gam_upper = [];
else
    gam_upper = Gam_NS*1.02; %16.655;%Gam_NS*(1.0);%4.389;
end

if isinf(gam_upper) 
    gam_upper = [];
end

% succeed = 1; % for indicating whether a feasible and stationary point has been found.
% tol = 5e-3*0.1; % for terminating the iteration when the change in value of optimization variables is small
% rou = 1e-4; % regulation parameter 
for hide = 1
if design_para.SwLogic == 0
    design_para.subset_weight = 1;
    design_para.full_search = 0;
    design_para.max_iter = 1;       
    design_para.ss_improve = 0;  
    theta_info.Theta1 = [theta_info.Theta1(1,1) theta_info.Theta1(end,end)];
    if gs_num == 2
        theta_info.Theta2 = [theta_info.Theta2(1,1) theta_info.Theta2(end,end)];
    end
end
F_theta = theta_fcn_info.F_theta;
d_F_theta = theta_fcn_info.d_F_theta;
FthetaNum = theta_fcn_info.FthetaNum;

FthetaNumX = FthetaNum; % used to indicate whether X is constant. FthetaNumX = 1 means X is constant
FthetaNumY = FthetaNum; 
switch design_para.XY_PD 
    case 1 %Y0 is constatnt, so define X          
        FthetaNumY = [1 0 0];
    case 2    %X0 is constatnt
        FthetaNumX = [1 0 0];      
    case 0
        FthetaNumX = [1 0 0]; 
        FthetaNumY = [1 0 0]; 
end
%% Plant dimension
n =  size(Gasym.A,1);
nw = size(Gasym.B1,2);
nu = size(Gasym.B2,2);
nz = size(Gasym.C1,1);
ny = size(Gasym.C2,1); 
nwz = n+nw+nz;
Inwz = eye(nwz);
In = eye(n);
%% regnum and regid determination,
regnum1 = size(theta_info.Theta1,1); %num. of subsets for theta1, only partition theta_1
regnum2 = size(theta_info.Theta2,1); % num. of subsets for theta2;
regnum = regnum1*regnum2; % Total num. of subsets. S shape from the bottom to order them, an example:
% 4 5 10
% 1 2 3
%% Determination of regid1 & regid2  
REGID = zeros(regnum,3);
for regid = 1:regnum
   if mod(regid,regnum1) == 0  
        regid1 = regnum1;
    else
        regid1 = mod(regid,regnum1);
   end
   if regid > regnum1
        regid2 = 1+ floor(regid/(regnum1+0.1));
    else
        regid2 = 1;
   end  
   REGID(regid,:) = [regid regid1 regid2]; 
end
%% Get constant part and linear part of plant matrices
Ga_00 = AugPltEval(Gasym,[0;0]);
A_0 = Ga_00.A; B1_0 = Ga_00.B1; B2_0 = Ga_00.B2;
C1_0 = Ga_00.C1; D11_0 = Ga_00.D11; D12_0 = Ga_00.D12;
C2_0 = Ga_00.C2; D21_0 = Ga_00.D21; D22_0 = Ga_00.D22;

Ga_10 = AugPltEval(Gasym,[1;0]);
Ga_m10 = AugPltEval(Gasym,[-1;0]);
A_1 = (Ga_10.A-Ga_m10.A)/2; B1_1 = (Ga_10.B1-Ga_m10.B1)/2; B2_1 = (Ga_10.B2-Ga_m10.B2)/2;
C1_1 = (Ga_10.C1-Ga_m10.C1)/2; D11_1 = (Ga_10.D11-Ga_m10.D11)/2; D12_1 = (Ga_10.D12-Ga_m10.D12)/2;
C2_1 = (Ga_10.C2-Ga_m10.C2)/2; D21_1 = (Ga_10.D21-Ga_m10.D21)/2; D22_1 = (Ga_10.D22-Ga_m10.D22)/2;

Ga_01 = AugPltEval(Gasym,[0;1]);
Ga_0m1 = AugPltEval(Gasym,[0;-1]);
A_2 = (Ga_01.A-Ga_0m1.A)/2; B1_2 = (Ga_01.B1-Ga_0m1.B1)/2; B2_2 = (Ga_01.B2-Ga_0m1.B2)/2;
C1_2 = (Ga_01.C1-Ga_0m1.C1)/2; D11_2 = (Ga_01.D11-Ga_0m1.D11)/2; D12_2 = (Ga_01.D12-Ga_0m1.D12)/2;
C2_2= (Ga_01.C2-Ga_0m1.C2)/2; D21_2 = (Ga_01.D21-Ga_0m1.D21)/2; D22_2 = (Ga_01.D22-Ga_0m1.D22)/2;


PD_order = 1;
% norm(M,'fro')> 1e-10 means M contains nonzero elements
if norm(Ga_10.A+Ga_m10.A-2*Ga_00.A,'fro')> 1e-10 || norm(Ga_01.A+Ga_0m1.A-2*Ga_00.A,'fro')> 1e-10 || norm(Ga_10.B1+Ga_m10.B1-2*Ga_00.B1,'fro') >1e-10 || norm(Ga_01.B1+Ga_0m1.B1-2*Ga_00.B1,'fro') >1e-10||...
        norm(Ga_10.C1+Ga_m10.C1-2*Ga_00.C1,'fro') >1e-10 || norm(Ga_01.C1+Ga_0m1.C1-2*Ga_00.C1,'fro') >1e-10
    PD_order = 2; % order of parameter dependence
    A_11 = (Ga_10.A+Ga_m10.A-2*Ga_00.A)/2; B1_11 = (Ga_10.B1+Ga_m10.B1-2*Ga_00.B1)/2; B2_11 = (Ga_10.B2+Ga_m10.B2-2*Ga_00.B2)/2;
    C1_11 = (Ga_10.C1+Ga_m10.C1-2*Ga_00.C1)/2; D11_11 = (Ga_10.D11+Ga_m10.D11-2*Ga_00.D11)/2; D12_11 = (Ga_10.D12+Ga_m10.D12-2*Ga_00.D12)/2;
    C2_11 = (Ga_10.C2+Ga_m10.C2-2*Ga_00.C2)/2; D21_11 = (Ga_10.D21+Ga_m10.D21-2*Ga_00.D21)/2; D22_11 = (Ga_10.D22+Ga_m10.D22-2*Ga_00.D22)/2;
    
    A_22 = (Ga_01.A+Ga_0m1.A-2*Ga_00.A)/2; B1_22 = (Ga_01.B1+Ga_0m1.B1-2*Ga_00.B1)/2; B2_22 = (Ga_01.B2+Ga_0m1.B2-2*Ga_00.B2)/2;
    C1_22 = (Ga_01.C1+Ga_0m1.C1-2*Ga_00.C1)/2; D11_22 = (Ga_01.D11+Ga_0m1.D11-2*Ga_00.D11)/2; D12_22 = (Ga_01.D12+Ga_0m1.D12-2*Ga_00.D12)/2;
    C2_22= (Ga_01.C2+Ga_0m1.C2-2*Ga_00.C2)/2; D21_22 = (Ga_01.D21+Ga_0m1.D21-2*Ga_00.D21)/2; D22_22 = (Ga_01.D22+Ga_0m1.D22-2*Ga_00.D22)/2;

    % get the cross_product term A_12, B1_12, C1_12
    Ga_11 = AugPltEval(Gasym,[1;1]);
    A_12 = Ga_11.A+Ga_00.A-Ga_10.A-Ga_01.A; 
    B1_12 = Ga_11.B1+Ga_00.B1-Ga_10.B1-Ga_01.B1; 
    C1_12 = Ga_11.C1+Ga_00.C1-Ga_10.C1-Ga_01.C1;     
end      
% Note that B2, C2, D12 and D21 are assumed to be affine in multi-convexity method
if PD_order == 1
    A_cell = {A_0,A_1,A_2}; 
    B1_cell = {B1_0,B1_1,B1_2; }; 
    B2_cell = {B2_0,B2_1,B2_2}; 
    C1_cell = {C1_0,C1_1,C1_2}; 
    D11_cell = {D11_0,D11_1,D11_2}; 
    D12_cell = {D12_0,D12_1,D12_2};
    C2_cell = {C2_0,C2_1,C2_2}; 
    D21_cell = {D21_0,D21_1,D21_2}; 
elseif PD_order == 2
    A_cell = {A_0,A_1,A_2;A_11,A_22,A_12}; % A_21 
    B1_cell = {B1_0,B1_1,B1_2; B1_11,B1_22,B1_12}; 
    B2_cell = {B2_0,B2_1,B2_2};  % B2 assumed to be affine
    C1_cell = {C1_0,C1_1,C1_2;C1_11,C1_22,C1_12}; 
    D11_cell = {D11_0,D11_1,D11_2;D11_11,D11_22,[]}; 
    D12_cell = {D12_0,D12_1,D12_2};% D12 assumed to be affine
    C2_cell = {C2_0,C2_1,C2_2}; % C2 assumed to be affine
    D21_cell = {D21_0,D21_1,D21_2}; % D21 assumed to be affine
end
    
Ga_ctmat.A_0 = A_0; Ga_ctmat.A_1 = A_1; Ga_ctmat.A_2 = A_2;
Ga_ctmat.B1_0 = B1_0; Ga_ctmat.B1_1 = B1_1; Ga_ctmat.B1_2 = B1_2;
Ga_ctmat.B2_0 = B2_0; Ga_ctmat.B2_1 = B2_1; Ga_ctmat.B2_2 = B2_2;
Ga_ctmat.C1_0 = C1_0; Ga_ctmat.C1_1 = C1_1; Ga_ctmat.C1_2 = C1_2;
Ga_ctmat.C2_0 = C2_0; Ga_ctmat.C2_1 = C2_1; Ga_ctmat.C2_2 = C2_2;
Ga_ctmat.D11_0 = D11_0; Ga_ctmat.D11_1 = D11_1; Ga_ctmat.D11_2 = D11_2; 
Ga_ctmat.D12_0 = D12_0; Ga_ctmat.D12_1 = D12_1; Ga_ctmat.D12_2 = D12_2;
Ga_ctmat.D21_0 = D21_0; Ga_ctmat.D21_1 = D21_1; Ga_ctmat.D21_2 = D21_2; 
Ga_ctmat.D22_0 = D22_0; Ga_ctmat.D22_1 = D22_1; Ga_ctmat.D22_2 = D22_2;

Ga_ctmat.A_cell = A_cell;
Ga_ctmat.B1_cell = B1_cell; Ga_ctmat.B2_cell = B2_cell;
Ga_ctmat.C1_cell = C1_cell; Ga_ctmat.C2_cell = C2_cell;
Ga_ctmat.D11_cell = D11_cell; Ga_ctmat.D12_cell = D12_cell;
Ga_ctmat.D21_cell = D21_cell; 

design_para.REGID = REGID;
design_para.gs_num = gs_num;
design_para.FthetaNum = FthetaNum;
design_para.FthetaNumX = FthetaNumX;
design_para.FthetaNumY = FthetaNumY;
design_para.F_theta= F_theta;
design_para.d_F_theta= d_F_theta;
design_para.gam_upper = gam_upper;
end