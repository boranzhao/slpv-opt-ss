function [opt,sol] = SLPV_Proj_Affine_MC(Gasym,Ga_ctmat,theta_info, design_para,opt_lmilab,opt_k)
%% Switching LPV controller design for LPV plants with affine parameter-dependence 
% Multi-convexity concept [Apkarian 2000] is used to solve the parameter-dependent LMIs  
% The default GS parameter number is 2, but can also be used for 1 GS parameter, in this case, the following settings should be imposed:
%   Theta2T = {0}, Delta2T = {0}, d_Thetah ={XX,0}, and Fcn_theta is only function of theta1, for instance Fcn_theta = @(x) [1 x(1)]  
% Modified on 10/07/2015 to use same constant X or Y, so that X(i) = X(j) need not be imposed. In this case, under average dwell time switching,
% constant terms are also imposed to be equal at switching surfaces. The contents for ASWCond_Ceq = 0 does not work now. 
%% Revision history
% 02/04/2017 ver 4: Revise the code especially the MC constraints to handle X Y containing second order
% term; add the penalty for overlapped regions between adjacent subsets
% 01/19/2017 ver 3: use different sig for main LMIs X,Y even for same
% subset
% 12/12/2016 ver2:
% use results from previous iterations as the initial value for
% current iteration for iteration on switching surface conditions. 
% Revised to deal with plants with second-order polynomial
% parameter dependence using Lyapunov variables with affine parameter
% dependence; furthermore, B2, C2, D12 and D21 are assumed to be affine
% w.r.t theta
% 12/07/2016 ver1: Revised to include the option of whether to use lambda & mu
% 11/10/2015 ver0: Created  

% Gasym, generalized plant, created using Matlab symbols
% theta_info a struct includes the following element:
%       ThetaT: cell array, gridded points of theta1 for all regions, ThetaT{1}  for subregion 1, Theta{2} for subregion 2,...
%       Theta2T: cell arry,gridded points of theta2 for all regions
%       Theta1: the variation range of theta1
%       Theta2: the variation range of theta2
%       d_Theta: cell array, bounds for derivatives of theta, d_Theta{1} for theta1, d_Theta{2} for theta2
% design_para a struct includes the following element:
%       REGID: [regid, regid1, regid2]
%       SwLogic: 0 for non-switching, 1 for hysteresis switching
%       XY_PD: 0 for constant X and Y, 1 for PDX and constant Y, 2 for constant X
%           and PD Y, 3 for PD X and Y
%       FthetaNum: a vector, each element show the number of constant matrices and matrices as a function of GS parameters theta1, theta2, ... 
%            for instance, [1 1 1] means one constant matrix, one matrix as a function of theta1, one matrix as a function of theta2,
%            while [1 2] means one constant matrix, two matrices as a function of theta1, no theta2.
%       Fcn_theta: a function handle, Ftheta(theta) will give all the scalar functions for
%                   the matrix variables, for instance in X = X0+
%                   f1(theta)X1+f2(theta)X2, Fcn_theta(theta) = [1 f1(theta) f2(theta)]
%       d_Fcn_theta: a function handle for the derivative of the scalar functions
% opt_yalmip: settings for yalmip

%% Paras 
if nargin < 6 || isempty(opt_k)
    X_k = []; Y_k = [];
    opt_k = [];
else 
    X_k = opt_k.X; Y_k = opt_k.Y;
end
xinit = [];
Theta1 = theta_info.Theta1;
Theta2 = theta_info.Theta2;
d_Theta = theta_info.d_Theta;
gs_num = theta_info.gs_num;
REGID = design_para.REGID;
SwLogic = design_para.SwLogic;
XY_PD = design_para.XY_PD;
FthetaNum = design_para.FthetaNum;
FthetaNumX = design_para.FthetaNumX;
FthetaNumY = design_para.FthetaNumY;
F_theta= design_para.F_theta;
LMI_relax = design_para.LMI_relax(1); % use to relax LMIs to guarantee strict feasiblity and avoid numerical problems. 
gam_upper = design_para.gam_upper; % upper limit of Gam value for each subset, equal to gamma value given non-switching LPV controller 
gam_min_method = design_para.gam_min_method; 
sig_limit = design_para.sig_limit;
XY_norm_limit = design_para.XY_norm_limit;
num_cost_gain = design_para.num_cost_gain;
% mc_relax = design_para.mc_relax;
% eps_ri = 1e-3; % to guarantee a solution that is in the interior of the solution set 
A_cell = Ga_ctmat.A_cell;
B1_cell = Ga_ctmat.B1_cell; B2_cell = Ga_ctmat.B2_cell;
C1_cell = Ga_ctmat.C1_cell; C2_cell = Ga_ctmat.C2_cell;
D11_cell = Ga_ctmat.D11_cell; D12_cell = Ga_ctmat.D12_cell;
D21_cell = Ga_ctmat.D21_cell; 

if norm(Theta2) == 0 || gs_num == 1 %% only one GS para theta1
   gs_num = 1; % Gain scheduling parameter number
   d_Theta{2} = 0;
   if max(REGID(:,3)) > 1
       disp('Error!REGID is not for one GS parameter!');
       return;
   end   
else
    gs_num = 2;
end
if XY_PD == 0
    d_Theta = {0,0};
    disp('Derivative of Thetah is set to 0 due to use of constant Lyapunov function');
end
%% Generalized plant parameter
n =  size(Gasym.A,1);
nw = size(Gasym.B1,2); %nu = size(Gasym.B2,2);
nz = size(Gasym.C1,1); %ny = size(Gasym.C2,1); 
%In = eye(n);  Inwz = eye(nwz); 
%% subregion parameter 
regnum1 = max(REGID(:,2));
regnum2 = max(REGID(:,3));
regnum = max(REGID(:,1)); % Total num. of subsets. S shape from the bottom to order them, an example:
% 4 5 6
% 1 2 3
if gam_min_method == 2
    regid1 = REGID(1:regnum,2); regid2 = REGID(1:regnum,3);
    if gs_num == 1
        subset_size = abs(Theta1(regid1,2)-Theta1(regid1,1));
        total_size = Theta1(end,end)-Theta1(1,1);
    elseif gs_num == 2
        subset_size = abs(Theta1(regid1,2)-Theta1(regid1,1)).*abs(Theta2(regid2,2)-Theta2(regid2,1));
        total_size = (Theta1(end,end)-Theta1(1,1))*(Theta2(end,end)-Theta2(1,1));
    end
end
%% optimization problem definition 
% optimization variables
sig = zeros(2,regnum); gam = zeros(1,regnum); 
lam = zeros(gs_num+1,regnum); mu = zeros(gs_num+1,regnum);
setlmis([]);
Gam = lmivar(1,[1 1]); 
gam_upper_slack = lmivar(1,[1 1]); 
if ~isempty(opt_k)
    xinit = opt_k.Gam+0.001; 
end
for regid = 1:regnum
    for id_Ftheta=1:sum(FthetaNum)%
        switch XY_PD 
            case 1 %Y0 is constatnt, so define X
                X(id_Ftheta,regid) =lmivar(1,[n 1]);                
            case 2    %X0 is constatnt
                Y(id_Ftheta,regid)=lmivar(1,[n 1]);   
            case 3
                X(id_Ftheta,regid)=lmivar(1,[n 1]); 
                Y(id_Ftheta,regid)=lmivar(1,[n 1]);
                if ~isempty(opt_k)
                    xinit =[xinit; ma2ve(opt_k.X(:,:,id_Ftheta,regid),1); ma2ve(opt_k.Y(:,:,id_Ftheta,regid),1)]; 
                end
        end
    end 
    gam(regid) = lmivar(1,[1 1]);
    sig(1,regid) = lmivar(1,[1 1]); 
    sig(2,regid) = sig(1,regid);    % set sig2 = sig1; 
%     sig(2,regid) =lmivar(1,[1 1]); %  sig(1,regid);%
    for i = 1: gs_num+1
        lam(i,regid) = lmivar(1,[1 1]);
        mu(i,regid) = lmivar(1,[1 1]);
        if ~(all(FthetaNumX<2) && all(FthetaNumY<2))
            eta(i,regid) = lmivar(1,[1 1]);
        end
    end 
    if ~isempty(opt_k) 
        if all(FthetaNumX<2) && all(FthetaNumY<2)
            tmp = [opt_k.lam(:,regid)'; opt_k.mu(:,regid)'];
        else
            tmp = [opt_k.lam(:,regid)'; opt_k.mu(:,regid)'; opt_k.eta(:,regid)'];
        end
        xinit =[xinit;opt_k.gam(regid);opt_k.sig(:,regid); tmp(:)];
    end
end      
switch XY_PD 
    case 1
        [Y(1,1),nt,sY] = lmivar(1,[n 1]);
        for regid = 2:regnum 
            Y(1,regid) = Y(1,1);% lmivar(3,sY);
        end
    case 2
        [X(1,1),nt,sX] = lmivar(1,[n 1]);
        for regid = 2:regnum 
            X(1,regid) = X(1,1); % lmivar(3,sX);
        end
    case 0
        [X(1,1),nt,sX]=lmivar(1,[n 1]); 
        [Y(1,1),nt,sY]=lmivar(1,[n 1]);         
        for regid = 2:regnum 
            X(1,regid) = X(1,1);% lmivar(3,sX);
            Y(1,regid) = Y(1,1); % lmivar(3,sY);
        end
        FthetaNumX = [1 0 0]; FthetaNumY = [1 0 0];
end
%% just for debugging
% for regid = 2: regnum
% %     X(:,regid) = X(:,1);
% %     Y(:,regid) = Y(:,1);
% %     sig(:,regid) = sig(:,1);
% %     lam(:,regid) = lam(:,1);
% %     mu(:,regid) = mu(:,1);
% end

if design_para.ss_cvx.method == 1 && XY_PD == 3
    if design_para.ss_cond == 1 && isempty(Y_k)
        for regid = 2: regnum
            Y(:,regid) = Y(:,1);
        end
%         Y(:,regid) = kron(Y(:,1),ones(regnum-1));
    elseif  design_para.ss_cond == 2 && isempty(X_k)
        for regid = 2: regnum
            X(:,regid) = X(:,1);
        end
    end
end
lminum = 0;
for regid = 1:regnum 
    if XY_PD == 0           
        lminum = lminum + 1;                
        lmiterm([lminum 1 1 X(1,1)],1,-1);
        lmiterm([lminum 1 1 0],LMI_relax); % for making positive definite instead of semi-definite
        lmiterm([lminum 2 1 0],-1);
        lmiterm([lminum 2 2 Y(1,1)],1,-1);
        lmiterm([lminum 2 2 0],LMI_relax);  
    end        
    %% Get the grid points for admissible region   
    regid1 = REGID(regid,2); regid2 = REGID(regid,3);
    %% Main LMIs
    for theta1 = Theta1(regid1,:) %check the vertices are enough due to use of multi-convexity concept
        for theta2 = Theta2(regid2,:); 
            theta = [theta1;theta2];
            Ftheta = F_theta([theta1;theta2]);  
            %% LMI3 for positivity of Lyapunov function. 
            if XY_PD ~= 0
                if all(FthetaNumX<2) && all(FthetaNumY<2)   % Both X and Y containt only affine terms
                    lminum = lminum +1;
                    for Id_Ftheta=1:sum(FthetaNumX)
                        lmiterm([lminum 1 1 X(Id_Ftheta,regid)],Ftheta(Id_Ftheta),-1);                                 
                    end 
                    for Id_Ftheta=1:sum(FthetaNumY)
                        lmiterm([lminum 2 2 Y(Id_Ftheta,regid)],Ftheta(Id_Ftheta),-1);                                 
                    end 
                    lmiterm([lminum 1 1 0],LMI_relax); % for making positive definite instead of semi-definite
                    lmiterm([lminum 2 1 0],-1);    
                    lmiterm([lminum 2 2 0],LMI_relax);
                else                                    % X or Y containt quadratic terms
                    lminum = lminum +1;
                    for Id_Ftheta=1:sum(FthetaNumX)
                        lmiterm([lminum 1 1 X(Id_Ftheta,regid)],Ftheta(Id_Ftheta),-1);                                 
                    end 
                    for Id_Ftheta=1:sum(FthetaNumY)
                        lmiterm([lminum 2 2 Y(Id_Ftheta,regid)],Ftheta(Id_Ftheta),-1);                                 
                    end 
                    lmiterm([lminum 2 1 0],-1);                     
                    for pos = 1:2
                        lmiterm([lminum pos pos 0],LMI_relax); % for making positive definite instead of semi-definite
                        lmiterm([lminum pos pos eta(1,regid)], 1,1);
                        for i = 1:gs_num
                            lmiterm([lminum pos pos eta(i+1,regid)],theta(i)^2,1);
                        end    
                    end
                    % multi-convexity constraints if X or Y including
                    % second-order terms                      
                    for i=1:gs_num
                        lminum = lminum +1;
                        lmiterm([lminum 1 1 X(sum(FthetaNumX(1:1+i)),regid)],1,1);
                        lmiterm([lminum 2 2 Y(sum(FthetaNumY(1:1+i)),regid)],1,1);
                        lmiterm([lminum 1 1 eta(i+1,regid)],-1,1);
                        lmiterm([lminum 2 2 eta(i+1,regid)],-1,1);
                    end
                end        
            end                       
            Ga = AugPltEval(Gasym, theta);
            A = Ga.A;
            B1 = Ga.B1;B2 = Ga.B2;
            C1 = Ga.C1;C2 = Ga.C2;
            D11 = Ga.D11; D12 = Ga.D12;
            D21 = Ga.D21; 
            
            if XY_PD == 1
                % LMI2 : no need to check the vertices of d_theta 
                lminum = lminum+1;    
                for Id_Ftheta=1:sum(FthetaNumY)
                    lmiterm([lminum 1 1 Y(Id_Ftheta,regid)],Ftheta(Id_Ftheta)*A,1,'s');
                    lmiterm([lminum 2 1 Y(Id_Ftheta,regid)],Ftheta(Id_Ftheta)*C1,1);                                        
                end                                       
                % No d_Y
                
                lmiterm([lminum 1 1 sig(2,regid)],-1,B2*B2');

                lmiterm([lminum 2 1 sig(2,regid)],-1,D12*B2');                    
                lmiterm([lminum 2 2 gam(regid)],-1,1);
                lmiterm([lminum 2 2 sig(2,regid)],-1,D12*D12');                       

                lmiterm([lminum 3 1 0], B1');                    
                lmiterm([lminum 3 2 0], D11');

%                     regid = tmp;
                lmiterm([lminum 3 3 gam(regid)],-1,1);                  
                
                
                for pos = 1:3
%                         lmiterm([lminum pos pos mu(1,regid)], 1,1);% no need
                    for i = 1:gs_num
                        lmiterm([lminum pos pos mu(i+1,regid)],theta(i)^2,1);
                    end   
                end
            elseif XY_PD == 2
                % LMI1: no need to check the vertices of d_theta 
                lminum = lminum+1;    
                for Id_Ftheta=1:sum(FthetaNumX)
                    lmiterm([lminum 1 1 X(Id_Ftheta,regid)],Ftheta(Id_Ftheta)*A',1,'s');
                    lmiterm([lminum 2 1 X(Id_Ftheta,regid)],Ftheta(Id_Ftheta)*B1',1);                                        
                end                
                % No d_X  
                lmiterm([lminum 1 1 sig(1,regid)],-1,C2'*C2);                 

                lmiterm([lminum 2 1 sig(1,regid)],-1,D21'*C2);                    
                lmiterm([lminum 2 2 gam(regid)],-1,1);
                lmiterm([lminum 2 2 sig(1,regid)],-1,D21'*D21);  

                lmiterm([lminum 3 1 0], C1);                    
                lmiterm([lminum 3 2 0], D11);
%                     regid = tmp;  
                lmiterm([lminum 3 3 gam(regid)],-1,1);     
%                     tmp = regid;
%                     regid = 1;
                for pos = 1:3
%                         lmiterm([lminum pos pos lam(1,regid)], 1,1); % no
%                         need to add this 
                    for i = 1:gs_num
                        lmiterm([lminum pos pos lam(i+1,regid)],theta(i)^2,1);
                    end    
                end     
            end                
            
            
            for d_theta1 = d_Theta{1} % d_Thetah = {0,0} is imposed when  XZ_PD ==0, so no problem for using for loop                
                for  d_theta2 = d_Theta{2}% d_Thetah{2} = 0, for one GS para case                               
                    % LMI1 
%                     tmp = regid;
%                     regid = 1;
                    if XY_PD == 1 || XY_PD == 3
                        lminum = lminum+1;    
                        for Id_Ftheta=1:sum(FthetaNumX)
                            lmiterm([lminum 1 1 X(Id_Ftheta,regid)],Ftheta(Id_Ftheta)*A',1,'s');
                            lmiterm([lminum 2 1 X(Id_Ftheta,regid)],Ftheta(Id_Ftheta)*B1',1);                                        
                        end                                       
                        % d_X
                        for Id_Ftheta = 2:1+FthetaNumX(2)
                            lmiterm([lminum 1 1 X(Id_Ftheta,regid)],d_theta1,1);
    %                                         -d_X = -(d_X+X_theta(:,:,Id_Ftheta,regid)*d_Fthetah(Id_Ftheta)*d_theta1); 
                        end                                 
                        for Id_Ftheta = 2+FthetaNumX(2):sum(FthetaNumX)%
                            lmiterm([lminum 1 1 X(Id_Ftheta,regid)],d_theta2,1);
                        end        
                        lmiterm([lminum 1 1 sig(1,regid)],-1,C2'*C2);                 

                        lmiterm([lminum 2 1 sig(1,regid)],-1,D21'*C2);                    
                        lmiterm([lminum 2 2 gam(regid)],-1,1);
                        lmiterm([lminum 2 2 sig(1,regid)],-1,D21'*D21);  

                        lmiterm([lminum 3 1 0], C1);                    
                        lmiterm([lminum 3 2 0], D11);
    %                     regid = tmp;  
                        lmiterm([lminum 3 3 gam(regid)],-1,1);                     

    %                     tmp = regid;
    %                     regid = 1;  
                        for pos = 1:3
    %                         lmiterm([lminum pos pos lam(1,regid)], 1,1); % no
    %                         need to add this 
                            for i = 1:gs_num
                                lmiterm([lminum pos pos lam(i+1,regid)],theta(i)^2,1);
                            end    
                        end
                    end
                    
                    if XY_PD == 2 || XY_PD == 3
                        % LMI2 
                        lminum = lminum+1;    
                        for Id_Ftheta=1:sum(FthetaNumY)
                            lmiterm([lminum 1 1 Y(Id_Ftheta,regid)],Ftheta(Id_Ftheta)*A,1,'s');
                            lmiterm([lminum 2 1 Y(Id_Ftheta,regid)],Ftheta(Id_Ftheta)*C1,1);                                        
                        end                                       
                        % d_Y
                        for Id_Ftheta = 2:1+FthetaNumY(2)
                            lmiterm([lminum 1 1 Y(Id_Ftheta,regid)],d_theta1,-1);
    %                                         -d_X = -(d_X+X_theta(:,:,Id_Ftheta,regid)*d_Fthetah(Id_Ftheta)*d_theta1); 
                        end                                 
                        for Id_Ftheta = 2+FthetaNumY(2):sum(FthetaNumY)%
                            lmiterm([lminum 1 1 Y(Id_Ftheta,regid)],d_theta2,-1);
                        end     
                        lmiterm([lminum 1 1 sig(2,regid)],-1,B2*B2');

                        lmiterm([lminum 2 1 sig(2,regid)],-1,D12*B2');                    
                        lmiterm([lminum 2 2 gam(regid)],-1,1);
                        lmiterm([lminum 2 2 sig(2,regid)],-1,D12*D12');                       

                        lmiterm([lminum 3 1 0], B1');                    
                        lmiterm([lminum 3 2 0], D11');

    %                     regid = tmp;
                        lmiterm([lminum 3 3 gam(regid)],-1,1);  
                        
                        for pos = 1:3
    %                         lmiterm([lminum pos pos mu(1,regid)], 1,1);% no need
                            for i = 1:gs_num
                                lmiterm([lminum pos pos mu(i+1,regid)],theta(i)^2,1);
                            end   
                        end
                    end
%                     if mc_relax == 1                        
%                     end                    
                end % d_theta2
            end % d_theta1          
        end %theta2
    end %theta1 
    
    % nonegativity of scalars
    if design_para.mc_relax == 1
        for i=1:gs_num+1
            lminum = lminum + 1;
            lmiterm([lminum 1 1 lam(i,regid)], -1,1);
            lminum = lminum + 1;
            lmiterm([lminum 1 1 mu(i,regid)], -1,1);  
            if ~(all(FthetaNumX<2) && all(FthetaNumY<2))
                lminum = lminum + 1;
                lmiterm([lminum 1 1 eta(i,regid)], -1,1);  
            end
        end
    end

    %% imposing convexity in direction theta_i 
    for hide = 1:1
    Theta_total = {Theta1(regid1,:); Theta2(regid2,:)};
    if  ~(all(FthetaNumX<2) && all(FthetaNumY<2)) % either X or Y contain quadratic terms
        for i = 1:gs_num 
            if (size(A_cell,1) < 2 || norm(A_cell{2,i},'fro')<1e-10 && norm(B1_cell{2,i},'fro')<1e-10)  && FthetaNumX(1+i)<2% no need to check every theta_i, if both A_cell{2,i} and B1_cell{2,i} are zero matrix
                theta_i_set = 0;
            else
                theta_i_set = Theta_total{i}; 
            end
            for theta_i = theta_i_set
                if gs_num == 1
                    j_index = 0;
                else
                    j_index = [1:i-1 i+1:gs_num];
                end
                for j = j_index
                    if gs_num == 1 || ((size(A_cell,1) < 2 || norm(A_cell{2,i},'fro')<1e-10 && norm(B1_cell{2,i},'fro')<1e-10) &&  FthetaNumX(1+i)<2) % no need to check every theta_j, if both A_cell{2,i} and B1_cell{2,i} are zero matrix
                        theta_j_set = 0;
                    else
                        theta_j_set = Theta_total{j};
                    end
                    for theta_j = theta_j_set
                        lminum = lminum + 1;
                        lmiterm([-lminum 1 1 X(sum(FthetaNumX(1:1+i)),regid)], 1, A_cell{1,1},'s'); % The 2nd row of A_cell is for quare terms 
                        lmiterm([-lminum 1 1 X(sum(FthetaNumX(1:i))+1,regid)], 1, A_cell{1,i+1},'s');   
%                         lmiterm([-lminum 1 1 X(sum(FthetaNumX(1:1+i)),regid)],theta_j, A_cell{1,j+1},'s'); 
                        lmiterm([-lminum 1 1 X(sum(FthetaNumX(1:1+i)),regid)],3*theta_i, A_cell{1,i+1},'s'); 
                        
                        if size(A_cell,1) >= 2
                            lmiterm([-lminum 1 1 X(1,regid)], 1, A_cell{2,i},'s'); % The 2nd row of A_cell is for quare terms                         
%                         lmiterm([-lminum 1 1 X(sum(FthetaNumX(1:j))+1,regid)], theta_j, A_cell{2,i},'s');                        
                            lmiterm([-lminum 1 1 X(sum(FthetaNumX(1:i))+1,regid)], 3*theta_i, A_cell{2,i},'s'); 
                          
%                             lmiterm([-lminum 1 1 X(sum(FthetaNumX(1:1+i)),regid)], theta_j^2, A_cell{2,j},'s');
%                         lmiterm([-lminum 1 1 X(sum(FthetaNumX(1:j))+1,regid)], theta_j^2, A_cell{2,i},'s'); 
                            lmiterm([-lminum 1 1 X(sum(FthetaNumX(1:1+i)),regid)], 6*theta_i^2, A_cell{2,i},'s'); 
                        end
                        
%                         
                        if size(A_cell,1)>=2 && norm(A_cell{2,3},'fro')> 1e-10 % The (2,3) element of A_cell is for cross-product terms
                            disp('A matrix should not contain cross-product term');
                        end
                        lmiterm([-lminum 1 1 sig(1,regid)],-C2_cell{i+1}'*C2_cell{i+1},1);

                        lmiterm([-lminum 2 1 X(sum(FthetaNumX(1:1+i)),regid)], B1_cell{1,1}',1); % The 2nd row of A_cell is for quare terms 
                        lmiterm([-lminum 2 1 X(sum(FthetaNumX(1:i))+1,regid)],B1_cell{1,i+1}', 1); 
%                         lmiterm([-lminum 2 1 X(sum(FthetaNumX(1:1+i)),regid)],B1_cell{1,j+1}',theta_j); 
                        lmiterm([-lminum 2 1 X(sum(FthetaNumX(1:1+i)),regid)], B1_cell{1,i+1}',3*theta_i); 
                        
                        if size(A_cell,1) >= 2
                            lmiterm([-lminum 2 1 X(1,regid)], B1_cell{2,i}', 1); % The 2nd row of B1_cell is for quare terms                         
    %                         lmiterm([-lminum 2 1 X(sum(FthetaNumX(1:j))+1,regid)], B1_cell{2,i}', theta_j);
                            lmiterm([-lminum 2 1 X(sum(FthetaNumX(1:i))+1,regid)], B1_cell{2,i}', 3*theta_i);  

    %                         lmiterm([-lminum 2 1 X(sum(FthetaNumX(1:1+i)),regid)], B1_cell{2,j}', theta_j^2);
    %                         lmiterm([-lminum 2 1 X(sum(FthetaNumX(1:j))+1,regid)], B1_cell{2,i}', theta_j^2); 
                            lmiterm([-lminum 2 1 X(sum(FthetaNumX(1:1+i)),regid)], B1_cell{2,i}', 6*theta_i^2); 
                        end
                        
                        if j~=0
                            lmiterm([-lminum 1 1 X(sum(FthetaNumX(1:1+i)),regid)],theta_j, A_cell{1,j+1},'s');                             
                            if size(A_cell,1) >= 2
                                lmiterm([-lminum 1 1 X(sum(FthetaNumX(1:j))+1,regid)], theta_j, A_cell{2,i},'s');
                                lmiterm([-lminum 1 1 X(sum(FthetaNumX(1:1+i)),regid)], theta_j^2, A_cell{2,j},'s');
                                lmiterm([-lminum 1 1 X(sum(FthetaNumX(1:j))+1,regid)], theta_j^2, A_cell{2,i},'s'); 
                            end
                            
                            lmiterm([-lminum 2 1 X(sum(FthetaNumX(1:1+i)),regid)],B1_cell{1,j+1}',theta_j); 
                            if size(B1_cell,1) >= 2
                                lmiterm([-lminum 2 1 X(sum(FthetaNumX(1:j))+1,regid)], B1_cell{2,i}', theta_j);
                                lmiterm([-lminum 2 1 X(sum(FthetaNumX(1:1+i)),regid)], B1_cell{2,j}', theta_j^2);
                                lmiterm([-lminum 2 1 X(sum(FthetaNumX(1:j))+1,regid)], B1_cell{2,i}', theta_j^2); 
                            end
                        end
                        
                        if size(B1_cell,1)>= 2 && norm(B1_cell{2,3},'fro')> 1e-10 
                            disp('B1 matrix should not contain cross-product term');
                        end            
                        lmiterm([-lminum 2 1 sig(1,regid)], -D21_cell{i+1}'*C2_cell{i+1},1);             
                        lmiterm([-lminum 2 2 sig(1,regid)], -D21_cell{i+1}'*D21_cell{i+1},1);
                    
                        lmiterm([-lminum 1 1 lam(i+1,regid)],1,1); 
                        lmiterm([-lminum 2 2 lam(i+1,regid)],1,1);    
                        if (size(C1_cell,1)>=2 && norm(C1_cell{2,i},'fro')> 1e-10)  || (size(D11_cell,1)>= 2 && norm(D11_cell{2,i},'fro')> 1e-10)
                            lmiterm([-lminum 3 1 0], C1_cell{2,i});
                            lmiterm([-lminum 3 2 0], D11_cell{2,i});
                            lmiterm([-lminum 3 3 lam(i+1,regid)],1,1);  
                        end
                    end
                end
            end
             if size(A_cell,1) >= 2
                % MC constraints from d^2f/dtheta_j^2
                if gs_num >= 2                
                    j_index = [1:i-1 i+1:gs_num];
                    for j = j_index
                        lminum = lminum +1;
                        lmiterm([lminum 1 1 X(sum(FthetaNumX(1:1+i)),regid)], 1, A_cell{2,j},'s'); 
                        lmiterm([lminum 1 1 X(sum(FthetaNumX(1:j))+1,regid)], 1, A_cell{2,i},'s');

                        lmiterm([lminum 2 1 X(sum(FthetaNumX(1:1+i)),regid)], B1_cell{2,j}',1);                
                        lmiterm([lminum 2 1 X(sum(FthetaNumX(1:j))+1,regid)], B1_cell{2,i}',1);               
                    end
                end
                % MC constraints from d^2f/dtheta_i^2
                lminum = lminum +1;
                lmiterm([lminum 1 1 X(sum(FthetaNumX(1:1+i)),regid)], 1, A_cell{2,i},'s');           
                lmiterm([lminum 2 1 X(sum(FthetaNumX(1:1+i)),regid)], B1_cell{2,i}',1);    
             end
        end
        
       
        for i = 1:gs_num 
            if (size(A_cell,1) < 2 || norm(A_cell{2,i},'fro')<1e-10 && norm(C1_cell{2,i},'fro')<1e-10) && FthetaNumY(1+i)<2 % no need to check every theta_i, if both A_cell{2,i} and C1_cell{2,i} are zero matrix
                theta_i_set = 0;
            else
                theta_i_set = Theta_total{i};
            end
            for theta_i = theta_i_set
                if gs_num == 1
                    j_index = 0;
                else
                    j_index = [1:i-1 i+1:gs_num];
                end
                for j = j_index
                    if gs_num == 1 || ((size(A_cell,1) < 2 || norm(A_cell{2,i},'fro')<1e-10 && norm(C1_cell{2,i},'fro')<1e-10) && FthetaNumY(1+i)<2) % no need to check every theta_j, if both A_cell{2,i} and C1_cell{2,i} are zero matrix
                        theta_j_set = 0;
                    else
                        theta_j_set = Theta_total{j};
                    end
                    for theta_j = theta_j_set
                        lminum = lminum + 1;
                        lmiterm([-lminum 1 1 Y(sum(FthetaNumY(1:1+i)),regid)],A_cell{1,1}, 1, 's'); % The 2nd row of A_cell is for quare terms 
                        lmiterm([-lminum 1 1 Y(sum(FthetaNumY(1:i))+1,regid)], A_cell{1,i+1}, 1,'s'); 
%                         lmiterm([-lminum 1 1 Y(sum(FthetaNumY(1:1+i)),regid)],A_cell{1,j+1},theta_j, 's'); 
                        lmiterm([-lminum 1 1 Y(sum(FthetaNumY(1:1+i)),regid)], A_cell{1,i+1},3*theta_i,'s'); 
                        if size(A_cell,1) >= 2
                            lmiterm([-lminum 1 1 Y(1,regid)], 1, A_cell{2,i},'s'); % The 2nd row of A_cell is for quare terms                         
    %                         lmiterm([-lminum 1 1 Y(sum(FthetaNumY(1:j))+1,regid)],A_cell{2,i}, theta_j, 's');
                            lmiterm([-lminum 1 1 Y(sum(FthetaNumY(1:i))+1,regid)], A_cell{2,i}, 3*theta_i,'s');  

    %                         lmiterm([-lminum 1 1 Y(sum(FthetaNumY(1:1+i)),regid)],A_cell{2,j}, theta_j^2, 's');
    %                         lmiterm([-lminum 1 1 Y(sum(FthetaNumY(1:j))+1,regid)],A_cell{2,i}, theta_j^2, 's'); 
                            lmiterm([-lminum 1 1 Y(sum(FthetaNumY(1:1+i)),regid)],A_cell{2,i}, 6*theta_i^2, 's'); 
                        end
                        
                        
                        lmiterm([-lminum 1 1 sig(2,regid)],-B2_cell{i+1}*B2_cell{i+1}',1);

                        lmiterm([-lminum 2 1 Y(sum(FthetaNumY(1:1+i)),regid)], C1_cell{1,1},1); % The 2nd row of A_cell is for quare terms 
                        lmiterm([-lminum 2 1 Y(sum(FthetaNumY(1:i))+1,regid)],C1_cell{1,i+1}, 1); 
%                         lmiterm([-lminum 2 1 Y(sum(FthetaNumY(1:1+i)),regid)],C1_cell{1,j+1},theta_j); 
                        lmiterm([-lminum 2 1 Y(sum(FthetaNumY(1:1+i)),regid)], C1_cell{1,i+1},3*theta_i); 
                        
                        if size(C1_cell,1) >= 2
                            lmiterm([-lminum 2 1 Y(1,regid)], C1_cell{2,i}, 1); % The 2nd row of C1_cell is for quare terms                         
    %                         lmiterm([-lminum 2 1 Y(sum(FthetaNumY(1:j))+1,regid)], C1_cell{2,i}, theta_j);
                            lmiterm([-lminum 2 1 Y(sum(FthetaNumY(1:i))+1,regid)], C1_cell{2,i}, 3*theta_i);  

    %                         lmiterm([-lminum 2 1 Y(sum(FthetaNumY(1:1+i)),regid)], C1_cell{2,j}, theta_j^2);
    %                         lmiterm([-lminum 2 1 Y(sum(FthetaNumY(1:j))+1,regid)], C1_cell{2,i}, theta_j^2); 
                            lmiterm([-lminum 2 1 Y(sum(FthetaNumY(1:1+i)),regid)], C1_cell{2,i}, 6*theta_i^2); 
                        end
                        
                        
                        if j~=0
                            lmiterm([-lminum 1 1 Y(sum(FthetaNumY(1:1+i)),regid)],A_cell{1,j+1},theta_j, 's');
                            if size(A_cell,1) >= 2
                                lmiterm([-lminum 1 1 Y(sum(FthetaNumY(1:j))+1,regid)],A_cell{2,i}, theta_j, 's');                            
                                lmiterm([-lminum 1 1 Y(sum(FthetaNumY(1:1+i)),regid)],A_cell{2,j}, theta_j^2, 's');
                                lmiterm([-lminum 1 1 Y(sum(FthetaNumY(1:j))+1,regid)],A_cell{2,i}, theta_j^2, 's'); 
                            end
                            
                            lmiterm([-lminum 2 1 Y(sum(FthetaNumY(1:1+i)),regid)],C1_cell{1,j+1},theta_j); 
                            if size(A_cell,1) >= 2
                                lmiterm([-lminum 2 1 Y(sum(FthetaNumY(1:j))+1,regid)], C1_cell{2,i}, theta_j);
                                lmiterm([-lminum 2 1 Y(sum(FthetaNumY(1:1+i)),regid)], C1_cell{2,j}, theta_j^2);
                                lmiterm([-lminum 2 1 Y(sum(FthetaNumY(1:j))+1,regid)], C1_cell{2,i}, theta_j^2);
                            end
                        end
                        
                        if size(C1_cell,1)>= 3 && norm(C1_cell{2,3},'fro')> 1e-10 
                            disp('C1 matrix should not contain cross-product term');                            
                        end            
                        lmiterm([-lminum 2 1 sig(2,regid)], -D12_cell{i+1}*B2_cell{i+1}',1);             
                        lmiterm([-lminum 2 2 sig(2,regid)], -D12_cell{i+1}*D12_cell{i+1}',1);
                        lmiterm([-lminum 1 1 mu(i+1,regid)],1,1); 
                        lmiterm([-lminum 2 2 mu(i+1,regid)],1,1); 
                        if (size(B1_cell,1)>=2 && norm(B1_cell{2,i},'fro')> 1e-10) || (size(D11_cell,1)>= 2 && norm(D11_cell{2,i},'fro')> 1e-10)
                            lmiterm([-lminum 3 1 0], B1_cell{2,i}');
                            lmiterm([-lminum 3 2 0], D11_cell{2,i}');
                            lmiterm([-lminum 3 3 mu(i+1,regid)],1,1);  
                        end
                    end
                end % j
            end % theta_i
            if size(A_cell,1) >= 2
                % MC constraints from d^2f/dtheta_j^2
                if gs_num >= 2
                    j_index = [1:i-1 i+1:gs_num];
                    for j = j_index
                        lminum = lminum +1;
                        lmiterm([lminum 1 1 Y(sum(FthetaNumY(1:1+i)),regid)],A_cell{2,j}, 1, 's'); 
                        lmiterm([lminum 1 1 Y(sum(FthetaNumY(1:j))+1,regid)], A_cell{2,i}, 1,'s');

                        lmiterm([lminum 2 1 Y(sum(FthetaNumY(1:1+i)),regid)], C1_cell{2,j}, 1);                
                        lmiterm([lminum 2 1 Y(sum(FthetaNumY(1:j))+1,regid)], C1_cell{2,i},1);               
                    end
                    % MC constraints from d^2f/dtheta_i^2
                    lminum = lminum +1;
                    lmiterm([lminum 1 1 Y(sum(FthetaNumY(1:1+i)),regid)], A_cell{2,i},1,'s');           
                    lmiterm([lminum 2 1 Y(sum(FthetaNumY(1:1+i)),regid)], B1_cell{2,i}',1);       
                end
            end
        end % i   
    elseif size(A_cell,1) >= 2 % Plant matrix A contains quadratic terms            
        for i = 1:gs_num 
            if (norm(A_cell{2,i},'fro')<1e-10 && norm(B1_cell{2,i},'fro')<1e-10) || XY_PD == 0 || XY_PD == 2% no need to check every theta_i, if both A_cell{2,i} and B1_cell{2,i} are zero matrix
                theta_i_set = 0;
            else
                theta_i_set = Theta_total{i}; 
            end
            for theta_i = theta_i_set
                if gs_num == 1
                    j_index = 0;
                else
                    j_index = [1:i-1 i+1:gs_num];
                end
                for j = j_index
                    if gs_num == 1 || (norm(A_cell{2,i},'fro')<1e-10 && norm(B1_cell{2,i},'fro')<1e-10) % no need to check every theta_j, if both A_cell{2,i} and B1_cell{2,i} are zero matrix
                        theta_j_set = 0;
                    else
                        theta_j_set = Theta_total{j};
                    end
                    for theta_j = theta_j_set
                        lminum = lminum + 1;
                        lmiterm([-lminum 1 1 X(1,regid)], 1, A_cell{2,i},'s'); % The 2nd row of A_cell is for quare terms 
                        if XY_PD == 1 || XY_PD == 3
                            lmiterm([-lminum 1 1 X(i+1,regid)], 1, A_cell{1,i+1},'s');  
                            lmiterm([-lminum 1 1 X(i+1,regid)], 3*theta_i, A_cell{2,i},'s');  
                            lmiterm([-lminum 1 1 X(j+1,regid)], theta_j, A_cell{2,i},'s');
                            if size(A_cell,1)>=2 && norm(A_cell{2,3},'fro')> 1e-10 % The (2,3) element of A_cell is for cross-product terms
                               lmiterm([-lminum 1 1 X(i+1,regid)], theta_j, A_cell{2,3},'s');
                            end
                        end
                        lmiterm([-lminum 1 1 sig(1,regid)],-C2_cell{i+1}'*C2_cell{i+1},1);

                        lmiterm([-lminum 2 1 X(1,regid)], B1_cell{2,i}',1);
                        if XY_PD == 1 || XY_PD == 3
                            lmiterm([-lminum 2 1 X(i+1,regid)], B1_cell{1,i+1}',1);
                            lmiterm([-lminum 2 1 X(i+1,regid)], B1_cell{2,i}',3*theta_i);
                            lmiterm([-lminum 2 1 X(j+1,regid)], B1_cell{2,i}',theta_j);
                            if size(B1_cell,1)>= 2 && norm(B1_cell{2,3},'fro')> 1e-10 
                                lmiterm([-lminum 1 1 X(i+1,regid)],B1_cell{2,3}',theta_j);
                            end
                        end
                        lmiterm([-lminum 2 1 sig(1,regid)], -D21_cell{i+1}'*C2_cell{i+1},1);             
                        lmiterm([-lminum 2 2 sig(1,regid)], -D21_cell{i+1}'*D21_cell{i+1},1);
                        lmiterm([-lminum 1 1 lam(i+1,regid)],1,1); 
                        lmiterm([-lminum 2 2 lam(i+1,regid)],1,1);    
                        if (size(C1_cell,1)>=2 && norm(C1_cell{2,i},'fro')> 1e-10)  || (size(D11_cell,1)>= 2 && norm(D11_cell{2,i},'fro')> 1e-10)
                            lmiterm([-lminum 3 1 0], C1_cell{2,i});
                            lmiterm([-lminum 3 2 0], D11_cell{2,i});
                            lmiterm([-lminum 3 3 lam(i+1,regid)],1,1);  
                        end
                    end
                end
            end
        end
        %
        for i = 1:gs_num 
            if (norm(A_cell{2,i},'fro')<1e-10 && norm(C1_cell{2,i},'fro')<1e-10) || XY_PD == 0 || XY_PD == 1  % no need to check every theta_i, if both A_cell{2,i} and C1_cell{2,i} are zero matrix
                theta_i_set = 0;
            else
                theta_i_set = Theta_total{i};
            end
            for theta_i = theta_i_set
                if gs_num == 1
                    j_index = 0;
                else
                    j_index = [1:i-1 i+1:gs_num];
                end
                for j = j_index
                    if gs_num == 1 || (norm(A_cell{2,i},'fro')<1e-10 && norm(C1_cell{2,i},'fro')<1e-10) % no need to check every theta_j, if both A_cell{2,i} and C1_cell{2,i} are zero matrix
                        theta_j_set = 0;
                    else
                        theta_j_set = Theta_total{j};
                    end
                    for theta_j = theta_j_set
                        lminum = lminum + 1;
                        lmiterm([-lminum 1 1 Y(1,regid)], 1, A_cell{2,i}','s'); % The 2nd row of A_cell is for quare terms                           
                        if XY_PD == 2 || XY_PD == 3
                            lmiterm([-lminum 1 1 Y(i+1,regid)], 1, A_cell{1,i+1}','s');                            
                            lmiterm([-lminum 1 1 Y(i+1,regid)], 3*theta_i, A_cell{2,i}','s');  
                            lmiterm([-lminum 1 1 Y(j+1,regid)], theta_j, A_cell{2,i}','s'); 
                            if size(A_cell,1)>=2 && norm(A_cell{2,3},'fro')> 1e-10   % The (2,3) of A_cell is for cross-product terms
                                lmiterm([-lminum 1 1 Y(i+1,regid)], theta_j, A_cell{2,3}','s');
                            end
                        end
                        lmiterm([-lminum 1 1 sig(2,regid)],-B2_cell{i+1}*B2_cell{i+1}',1);

                        lmiterm([-lminum 2 1 Y(1,regid)], C1_cell{2,i},1);
                        if XY_PD == 2 || XY_PD == 3
                            lmiterm([-lminum 2 1 Y(i+1,regid)], C1_cell{1,i+1},1);
                            lmiterm([-lminum 2 1 Y(i+1,regid)], C1_cell{2,i},3*theta_i);
                            lmiterm([-lminum 2 1 Y(j+1,regid)], C1_cell{2,i},theta_j);
                            if size(C1_cell,1)>= 3 && norm(C1_cell{2,3},'fro')> 1e-10 
                                lmiterm([-lminum 1 1 Y(i+1,regid)],C1_cell{2,3},theta_j);                            
                            end        
                        end
                        lmiterm([-lminum 2 1 sig(2,regid)], -D12_cell{i+1}*B2_cell{i+1}',1);             
                        lmiterm([-lminum 2 2 sig(2,regid)], -D12_cell{i+1}*D12_cell{i+1}',1);
                        lmiterm([-lminum 1 1 mu(i+1,regid)],1,1); 
                        lmiterm([-lminum 2 2 mu(i+1,regid)],1,1); 
                        if (size(B1_cell,1)>=2 && norm(B1_cell{2,i},'fro')> 1e-10) || (size(D11_cell,1)>= 2 && norm(D11_cell{2,i},'fro')> 1e-10)
                            lmiterm([-lminum 3 1 0], B1_cell{2,i}');
                            lmiterm([-lminum 3 2 0], D11_cell{2,i}');
                            lmiterm([-lminum 3 3 mu(i+1,regid)],1,1);  
                        end
                    end
                end % j
            end % theta_i
        end % i
    else % both plant matrices and X & Y contain only affine terms
        for i = 1:gs_num 
            % 1
            lminum = lminum + 1;
            if XY_PD == 1 || XY_PD == 3                
                lmiterm([-lminum 1 1 X(i+1,regid)], 1, A_cell{i+1},'s');  
                lmiterm([-lminum 2 1 X(i+1,regid)], B1_cell{i+1}',1);
            end
            lmiterm([-lminum 1 1 sig(1,regid)],-C2_cell{i+1}'*C2_cell{i+1},1);
            lmiterm([-lminum 2 1 sig(1,regid)], -D21_cell{i+1}'*C2_cell{i+1},1);  
            lmiterm([-lminum 2 2 sig(1,regid)], -D21_cell{i+1}'*D21_cell{i+1},1);
            lmiterm([-lminum 1 1 lam(i+1,regid)],1,1); 
            lmiterm([-lminum 2 2 lam(i+1,regid)],1,1); 
            
            % 2
             lminum = lminum + 1;
            if XY_PD == 2 || XY_PD == 3               
                lmiterm([-lminum 1 1 Y(i+1,regid)], 1, A_cell{i+1}','s');    
                lmiterm([-lminum 2 1 Y(i+1,regid)], C1_cell{i+1},1);                
            end
            lmiterm([-lminum 1 1 sig(2,regid)],-B2_cell{i+1}*B2_cell{i+1}',1);   
            lmiterm([-lminum 2 1 sig(2,regid)], -D12_cell{i+1}*B2_cell{i+1}',1); 
            lmiterm([-lminum 2 2 sig(2,regid)], -D12_cell{i+1}*D12_cell{i+1}',1);
            lmiterm([-lminum 1 1 mu(i+1,regid)],1,1);   
            lmiterm([-lminum 2 2 mu(i+1,regid)],1,1);
        end    
    end
    end
%  
    % cost function 
    if gam_min_method == 1
        lminum = lminum + 1;  
        lmiterm([lminum 1 1 gam(regid)],1,1);
        lmiterm([lminum 1 1 Gam],-1,1);
    elseif gam_min_method == 2 
        if ~isempty(gam_upper)
            lminum = lminum + 1;  
            lmiterm([lminum 1 1 gam(regid)],1,1);
            lmiterm([lminum 1 1 0],-gam_upper);
            lmiterm([lminum 1 1 gam_upper_slack],-1,1); % make the  constraint soft
        end
    end
end% regid

if gam_min_method == 2
    lminum = lminum + 1;
    for regid = 1: regnum
        if isempty(design_para.subset_weight)
            lmiterm([lminum 1 1 gam(regid)],subset_size(regid)*theta_info.agg_subset(regid)/sum(theta_info.agg_subset*subset_size),1); 
        else
            lmiterm([lminum 1 1 gam(regid)],design_para.subset_weight(regid),1);
        end
    end
    lmiterm([lminum 1 1 gam_upper_slack],100,1); % penalize violation of the constraints
    lmiterm([lminum 1 1 Gam],1,-1);
    
    lminum = lminum + 1;
    lmiterm([lminum 1 1 gam_upper_slack],-1,1);
end    
lmisys = getlmis;
%% LMIs for the monotonic property of Lyapunov function      
% practical validity constrainsts are imposed to convexify the
% Note that X Y follows the notations, different from the notations
% used in my Auotomatica Paper
% nonconvex switching surface conditions.
% switching surface index increases in the theta1 direction first
% in the way of 
% 4 5 6
% 1 2 3
% and then in theta2 direction in the way of 
% 2 4
% 1 3   
% Paras.XY_PD = XY_PD;
% Paras.eps_eq = eps_eq;  
% Paras.ss_cond = ss_cond;
%% for switching surfaces in theta1 direction
%  SwLogic = 0;
if SwLogic ~= 0  
for regid2 = 1:regnum2
    RegId1 = REGID(find(REGID(:,3)==regid2),2);
    RegId1 = RegId1'; RegId1 = sort(RegId1);     
    RegId = (regid2-1)*regnum1+RegId1;
    for id1 = RegId1(1:end-1)    
        for theta2 = Theta2(regid2,:)%becasue inequalities is affine w.r.t theta1 & theta2, so check of vertices are enough
            % P(RegId(i)) >= P(RegId(i+1)) @ theta2 = Theta2(i,2)
            theta1 = Theta1(id1,end);                                                      
            Ftheta = F_theta([theta1;theta2]);   
            jr = [RegId(id1) RegId(id1+1)];
            [lmisys,lminum] = SS_LMIs_CtlDsgn(lmisys,X,Y,Ftheta,jr,design_para,lminum,X_k,Y_k);

            % P(RegId(i+1)) >= P(RegId(i)) @ theta2 = Theta2(i+1,1)
            theta1 = Theta1(id1+1,1);
            Ftheta = F_theta([theta1;theta2]);                                   
            jr = [RegId(id1+1) RegId(id1)];                            
            [lmisys,lminum] = SS_LMIs_CtlDsgn(lmisys,X,Y,Ftheta,jr,design_para,lminum,X_k,Y_k);
        end 
    end           
end     
%% for switching surfaces in theta2 direction 
for regid1 = 1:regnum1
    RegId2 = REGID(find(REGID(:,2)==regid1),3);
    RegId2 = RegId2';RegId2 = sort(RegId2);       
    RegId = (RegId2-1)*regnum1+regid1;
    for id2 = RegId2(1:end-1) 
        for theta1 = Theta1(regid1,:) %becasue inequalities is affine w.r.t theta1 & theta2, so check of vertices are enough
            % P(RegId(i)) >= P(RegId(i+1)) @ theta2 = Theta2(i,2)
            theta2 = Theta2(id2,end);
            Ftheta = F_theta([theta1;theta2]);   
            jr = [RegId(id2) RegId(id2+1)];
            [lmisys,lminum] = SS_LMIs_CtlDsgn(lmisys,X,Y,Ftheta,jr,design_para,lminum,X_k,Y_k);            
            % P(RegId(i+1)) >= P(RegId(i)) @ theta2 = Theta2(i+1,1)
            theta2 = Theta2(id2+1,1);
            Ftheta = F_theta([theta1;theta2]);  
            jr = [RegId(id2+1) RegId(id2)];      
            [lmisys,lminum] = SS_LMIs_CtlDsgn(lmisys,X,Y,Ftheta,jr,design_para,lminum,X_k,Y_k);  
        end           
    end
end
end   % if SwLogic   

setlmis(lmisys);
if sum(num_cost_gain) ~=0        
    lminum_num_stab = lminum + 1; lminum = lminum_num_stab;       
    if num_cost_gain(1) ~= 0
        XY_norm_var = lmivar(1,[1 1]);
        if ~isempty(opt_k) 
            xinit =[xinit;(opt_k.obj-opt_k.Gam)/(num_cost_gain(1)+1e-10)]; 
        end
        for regid = 1:1 %? regnum        
%             % limit the 2-norm of X0^{j} Y0^{j}:limiting the eigenvalue
%             lminum = lminum+1;
%             lmiterm([lminum 1 1 X(1,regid)],1,1);
%             lmiterm([lminum 1 1 XY_norm_var],-1,1);
%             lmiterm([lminum 1 1 0],-XY_norm_limit);%
% 
%             lminum = lminum+1;
%             lmiterm([lminum 1 1 Y(1,regid)],1,1);
%             lmiterm([lminum 1 1 XY_norm_var],-1,1);
%             lmiterm([lminum 1 1 0],-XY_norm_limit);%   

            % limit the 2-norm of X0^{j} Y0^{j}         
            % X'*X <= (XY_norm_var + XY_norm_limit)^2 
            % [XY_norm_var + XY_norm_limit X'; X XY_norm_var + XY_norm_limit]
            % >= 0
            for Id_Ftheta = 1: sum(FthetaNumX) % for test 
                lminum = lminum+1;
                lmiterm([-lminum 1 1 XY_norm_var],1,1);
                lmiterm([-lminum 1 1 0],XY_norm_limit);%
                lmiterm([-lminum 2 1 X(Id_Ftheta,regid)],1,1);
                lmiterm([-lminum 2 2 XY_norm_var],1,1);
                lmiterm([-lminum 2 2 0],XY_norm_limit);%
            end
            
            for Id_Ftheta = 1: sum(FthetaNumY) % for test 
                lminum = lminum+1;
                lmiterm([-lminum 1 1 XY_norm_var],1,1);
                lmiterm([-lminum 1 1 0],XY_norm_limit);% 
                lmiterm([-lminum 2 1 Y(Id_Ftheta,regid)],1,1); 
                lmiterm([-lminum 2 2 XY_norm_var],1,1);
                lmiterm([-lminum 2 2 0],XY_norm_limit);%
            end
        end
        lmiterm([lminum_num_stab 1 1 XY_norm_var],num_cost_gain(1),1);
        
        lminum = lminum+1; %X_norm_var >= 0 
        lmiterm([lminum 1 1 XY_norm_var],-1,1);
    end

    if num_cost_gain(2) ~= 0
        sig_var = lmivar(1,[1 1]); 
        if ~isempty(opt_k) 
            xinit =[xinit;(opt_k.obj-opt_k.Gam)/(num_cost_gain(2)+1e-10)]; 
        end
        for regid = 1:regnum      
            lminum = lminum + 1;
            lmiterm([lminum 1 1 sig(1,regid)],1,1);
            lmiterm([lminum 1 1 sig_var],-1,1);
            lmiterm([lminum 1 1 0],-sig_limit);%     
            
            lminum = lminum + 1;
            lmiterm([lminum 1 1 sig(2,regid)],1,1);
            lmiterm([lminum 1 1 sig_var],-1,1);
            lmiterm([lminum 1 1 0],-sig_limit);%     
        end         
        lmiterm([lminum_num_stab 1 1 sig_var],num_cost_gain(2),1);
        lminum = lminum+1; %sig_var >= 0 
        lmiterm([lminum 1 1 sig_var],-1,1);
    end    
    obj = lmivar(1,[1 1]); 
    if ~isempty(opt_k) 
        xinit =[xinit;opt_k.obj+3*(opt_k.obj-opt_k.Gam)+0.002]; 
    end
    lmiterm([lminum_num_stab 1 1 obj],-1,1);
    lmiterm([lminum_num_stab 1 1 Gam],1,1);   
end
if design_para.overlap_gain ~= 0
    lmiterm([lminum_num_stab 1 1 0], max(sum(subset_size)-total_size-design_para.overlap_limit,0)*design_para.overlap_gain);   
end

lmisys = getlmis;        
nvar = decnbr(lmisys); %number of decision variable.
c = zeros(nvar,1);
if sum(num_cost_gain) == 0
    c(1)=1;
else
    c(end) = 1;
end 
if isempty(xinit)
    [fopt,xopt] = mincx(lmisys,c,opt_lmilab);
else
    [fopt,xopt] = mincx(lmisys,c,opt_lmilab,xinit);
end

if ~isempty(fopt)
%     %% check whether constraints are satisfied
%     lmi_succeed = true;
%     evals = evallmi(lmisys,xopt);        
%     for i = 1:lminum
%         [lhs,rhs] = showlmi(evals,i);
%         if ~all(real(eig(lhs-rhs))<= 0)               
%             lmi_succeed = false;max(real(eig(lhs-rhs)))
%     %                 pause;
%     %                 break;
%         end
%     end 
%     if lmi_succeed
%         disp('All constraints satisfied!');            
%     else
%         disp('Oops!Not all the constraints satisfied!');
%     end    
    %% Extract the value of optimization variables
    opt.X = zeros(n,n,sum(FthetaNumX),regnum);
    opt.Y = zeros(n,n,sum(FthetaNumY),regnum);
    opt.lam = zeros(gs_num+1,regnum);
    opt.mu = zeros(gs_num+1,regnum);
    opt.gam = zeros(1,regnum);
    opt.sig = zeros(2,regnum);
    for regid = 1:regnum
        for id_Ftheta = 1:sum(FthetaNumX) 
            opt.X(:,:,id_Ftheta,regid) = dec2mat(lmisys,xopt,X(id_Ftheta,regid));
        end
        for id_Ftheta = 1:sum(FthetaNumY) 
            opt.Y(:,:,id_Ftheta,regid) = dec2mat(lmisys,xopt,Y(id_Ftheta,regid));
        end
        for i = 1:gs_num+1
            opt.lam(i,regid) = dec2mat(lmisys,xopt,lam(i,regid));
            opt.mu(i,regid) = dec2mat(lmisys,xopt,mu(i,regid));
        end
        opt.gam(regid)  = dec2mat(lmisys,xopt,gam(regid));
        opt.sig(1,regid)  = dec2mat(lmisys,xopt,sig(1,regid));
        opt.sig(2,regid)  = dec2mat(lmisys,xopt,sig(2,regid));
    end   

    opt.Gam =  dec2mat(lmisys,xopt,Gam);
    if exist('obj','var') == 1
        opt.obj = dec2mat(lmisys,xopt,obj);
    else
        opt.obj = NaN;
    end
else
    opt.X = []; opt.Y = []; opt.gam = inf*ones(1,regnum); opt.sig = [];
    opt.lam =[]; opt.mu = []; opt.Gam = inf; opt.obj = inf;
end
 sol = [];    