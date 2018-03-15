function [Gasym,theta_fcn_info,gs_num,theta_info,num_gain] = plant_sel(plant_index,Theta1,gs_num)

%% plant_index selection
    
if plant_index == 1 % numerical example ?originally from Sato 2011?
    for hide = 1
    num_gain = 1; % for imporving numerical stability 
    syms theta1
    Gasym.A = [-4 3 5;0 7 -5;0.1 -2 -3]+theta1*[1 0 1;2 0 -5;2 5 1.5]/num_gain;%+ theta1^2*eye(3)*0.5/num_gain; 
    Gasym.B1 = [1;-2;1]; % +theta1^2*[0.2;-0.5;0.2];
    if gs_num == 1
        Gasym.B2 = [0;16;-10]+theta1*[1;-5;3.5]/num_gain;
    elseif gs_num == 2
        syms theta2;
        Gasym.B2 = [0;16;-10]+theta2*[1;-5;3.5]/num_gain;
    end
    Gasym.C1 = [1 1 0];
    Gasym.D11 = 0;
    Gasym.D12 = 1;
    Gasym.C2 = [0 1 0];
    Gasym.D21 = 2;
    Gasym.D22 = 0;
    
    %% poster filter u to remove the parameter dependence of B2
%     n =  size(Gasym.A,1);
%     nw = size(Gasym.B1,2);
%     nz = size(Gasym.C1,1);
%     ny = size(Gasym.C2,1);
%     Fu = ss(tf(1,[1/5e3 1])); n_fu = order(Fu);
%     Gasym1.A = [Gasym.A Gasym.B2*Fu.c; zeros(1,n) Fu.a];
%     Gasym1.B1 = [Gasym.B1; zeros(1,nw)];
%     Gasym1.B2 = [zeros(n,size(Fu.b,2)); Fu.b];
%     Gasym1.C1 = [Gasym.C1 Gasym.D12*Fu.c]; 
%     Gasym1.D11 = Gasym.D11; Gasym1.D12 = zeros(nz,size(Fu.b,2));
%     Gasym1.C2 = [Gasym.C2 zeros(ny,n_fu)];
%     Gasym1.D21 = Gasym.D21; Gasym1.D22 = Gasym.D22;     
%     Gasym = Gasym1;
%         
     Theta1 = Theta1*num_gain;    
     d_theta = 1*num_gain;
    if gs_num == 1
        Theta2 =0 ;
        d_Theta = {[-d_theta d_theta],0};
        FthetaNum = [1 1];     
        F_theta = @(x) [1 x(1,1)]; %function for PD matrices
        d_F_theta = @(x) [0 1];  %function 
    elseif gs_num == 2
        Theta2 =  [-1  -0.6351;-0.6451 1];%[-3 0.1;-0.1 3]*num_gain;
        d_Theta = {[-d_theta d_theta],[-d_theta d_theta]};
        FthetaNum = [1 1 1];     
        F_theta = @(x) [1 x(1) x(2)]; %function for PD matrices
        d_F_theta = @(x) [0 1 1];  %function
    end  
    ss_dist_min = 0.01;

%       %% Hinf controller design
%     Ga = AugPltEval(Gasym,0.5);
%     a = Ga.A;
% 
%     b = [Ga.B1 Ga.B2];
%     c = [Ga.C1;Ga.C2];
%     d = [Ga.D11 Ga.D12;
%         Ga.D21 Ga.D22];
%     ga = ss(a,b,c,d);    
%     [Kinf,CL,Gam,INFO] = hinfsyn(ga,1,1,'method','lmi')
%     return;    
    
    end    
elseif plant_index == 2 % engine control example
    %%
    num_gain = 1;
    XY_quadterm = 0;
    syms theta1 theta2 %  theta1: 1./(120./N+ c./ma);  theta2: engine speed;
%     theta1 = 3.0628;
%     theta2 = 2000;
    if gs_num == 1
        theta2 = 3000;
    end
    %% test the range of 1/T and its derivative.
%     N = 800:100:4000;
%     ma = 0.1:0.02:1;
%     N_dot = 2000;
%     ma_dot = 1;
%     c = 5.33e-2;
%     [N,ma] = meshgrid(N,ma);
%     Tinv =1./(120./N+ c./ma);
%     Tinv_dot = (120*N_dot./N.^2+c*ma_dot./ma.^2)./(120./N+c./ma).^2;
%     figure;surf(N,ma,Tinv);
%     figure;surf(N,ma,Tinv_dot);
%     Tinv_max = max(max(Tinv));
%     Tinv_min = min(min(Tinv));
%     
%     Tinv_dot_max = max(max(Tinv_dot));
%     Tinv_dot_min = min(min(Tinv_dot));

% some special value
% N = 4000, ma = 0.8; Tinv = 10.3493
% N = 1000, ma = 0.15: Tinv = 2.1038
    %% plant parameters      
    ncyl = 4;
    tau = 2*60*(ncyl-1)/ncyl/(theta2/num_gain);
    Rstoich = 14.7;
    g = Rstoich; 
%     Filter = ss(tf(1,[3e0 1]));
    
%     T = 1/theta1;
%     s = tf('s'); 
%     sysDelay = exp(-T*s);

%     Filter = ss(tf(1,[0.01 1]));
    delayMod = 1; % models for approximating pure delay: 
    % 0 for no delay
    % 1 for 1st order ; 
    % 2 for 2nd order app                                                                                                                      rox. from Marius; 
    % 3 for ordinary 2nd order 
        
    switch delayMod 
        case 0 % 1st order system without delay
            Ap = -1/tau;
            Bp = g;
            Cp = 1/tau;    
            Dp = 0;             
        case 1  % 1st order with delay     
            Ap = [-1/tau 4*g*theta1; 0 -2*theta1];
            Bp = [-g 1].';
            Cp = [1/tau 0];  
%             Bp = [-g 1]'/tau;
%             Cp = [1 0];   
            Dp =0;
        % higher order pade approx.
        %     plt = tf(g,[tau 1]);
        %     T = 1/theta1;   
        %     [num,den] = pade(T,3);
        %     delay = tf(num,den);
        %     [Ap,Bp,Cp,Dp] = ssdata(plt*delay);
        case 2 % 2nd order pade from Marius
            Ap = [-1/tau 6*theta1^2 -2*theta1; 0 0 1; 0 -6*theta1^2 -4*theta1];
%             Bp = [0 0 1/tau].';
%             Cp = [g 0 0];  
            Bp = [0 0 1].';
            Cp = [g 0 0]/tau; 
            Dp = 0;  
        case 3 % Ordinary 2nd order pade 
            Ap = [-1/tau 0 -12*theta1; 0 0 1; 0 -12*theta1^2 -6*theta1];
            Bp = [1 0 1].';
            Cp = [g/tau 0 0];     
            Dp = 0;
           %% to determine the frequency of the filter used to control a system with pure time delay 
        %     Tmax= 0.683060;
        %     [num,den] = pade(Tmax,1);
        %     delay_app = tf(num,den);
        %     s = tf('s');
        %     delay_real = exp(-Tmax*s);
        % %     delay_real = tf([-2*Tmax, 6],[Tmax^2,4*Tmax,6]);
        %     figure;
        %     h = bodeplot(delay_app,delay_real); % returns handle of the bode plot
        %     % set the following options for the bode plot
        %     setoptions(h,'PhaseMatching', 'on', 'PhaseMatchingValue', 0);
    end
%     Cp*inv(s*eye(length(Ap))-Ap)*Bp

      % new plant including the filter 
      if exist('Filter','var')
          n=length(Ap);
          plt.a =  [Filter.a Filter.b*Cp; zeros(n,1) Ap];
          plt.b = [Filter.b*Dp; Bp];
          plt.c = [Filter.c Filter.d*Cp];
          plt.d = Filter.d*Dp;
          Ap = plt.a; Bp = plt.b;Cp=plt.c; Dp = plt.d;
      end

      %%
    Wi = ss(tf([0.6 2],[1 0])); % 0.6,a first order system instead of an integrator to enforce D21 ~= 0
    assignin('base','Wi',Wi);
    %     Plt = ss(tf(1,[1 1]));   
%     Ap = Plt.a;
%     Bp = Plt.b;
%     Cp = Plt.c;
%     Dp = Plt.d;
    
%     % without delay
%     Ap = -1/tau;
%     Bp = g;
%     Cp = 1/tau;    
%     Dp = 0; 
%     plt = ss(Ap,Bp,Cp,Dp);
    np = size(Ap,1);
%     assignin('base','Plt',Plt);      
%% weigthting function 
    We = 1;nwe = 0;
    Wu = 1e1; nwu = 0;
    assignin('base','We',We);
    assignin('base','Wu',Wu);
    nint = 1; % order of integrator  
   
    Gasym.A = [Ap, zeros(np,nwe+nwu+nint);...           
        -Wi.b*Cp, zeros(nint,nwe+nwu) Wi.a];
    Gasym.B1 =  [zeros(np,2);...
        -Wi.b, Wi.b];
    Gasym.B2 = [Bp;-Wi.b*Dp];
    Gasym.C1 = [[-Wi.d*Cp zeros(1,nwe+nwu), Wi.c]*We;...
        zeros(1,np+nwe), zeros(1,nint)];
    Gasym.D11 = [-We*Wi.d, We*Wi.d;...
        zeros(1,2)];
    Gasym.D12 = [-We*Wi.d*Dp; Wu];
    Gasym.C2 = [-Wi.d*Cp zeros(1,nwe+nwu), Wi.c];  
    Gasym.D21 = [-Wi.d Wi.d];
    Gasym.D22 = -Wi.d*Dp;
    
    %% poster filter u to remove the parameter dependence of B2
%     n =  size(Gasym.A,1);
%     nw = size(Gasym.B1,2);
%     nz = size(Gasym.C1,1);
%     ny = size(Gasym.C2,1);
%     Fu = ss(tf(1,[1/5e3 1])); n_fu = order(Fu);
%     Gasym1.A = [Gasym.A Gasym.B2*Fu.c; zeros(1,n) Fu.a];
%     Gasym1.B1 = [Gasym.B1; zeros(1,nw)];
%     Gasym1.B2 = [zeros(n,size(Fu.b,2)); Fu.b];
%     Gasym1.C1 = [Gasym.C1 Gasym.D12*Fu.c]; 
%     Gasym1.D11 = Gasym.D11; Gasym1.D12 = zeros(nz,size(Fu.b,2));
%     Gasym1.C2 = [Gasym.C2 zeros(ny,n_fu)];
%     Gasym1.D21 = Gasym.D21; Gasym1.D22 = Gasym.D22;     
%     Gasym = Gasym1;
    
    if gs_num == 2
        Theta2 = [800  6000]*num_gain; % heuristic
%         Theta2 = [800 3500; 3000  6000]; % heuristic
%         Theta2 = [800 4403.4; 4390.3  6000];% pso_4subs_XY3_Proj_2perc, [800 3500; 3000  6000]/num_gain; % [800 3500; 3000 6000]
%         Theta2 = [800 4426.8; 2543.9  6000];
        d_Theta = {[-17.46 17.46],[-2000 2000]*num_gain}; %Droumax and Droumin  
%         d_Theta = {0,0/num_gain}; %Droumax and Droumin  
        if XY_quadterm == 0 % without second order term
            FthetaNum = [1 1 1];
            F_theta = @(x) [1 x(1) x(2)]; %function for PD matrices
            d_F_theta = @(x) [0 1 1];  %function for derivative of PD matrices
        elseif XY_quadterm == 1
            FthetaNum = [1 2 1];
            F_theta = @(x) [1 x(1)  x(1)^2 x(2)]; %function for PD matrices
            d_F_theta = @(x) [0 1 2*x(1) 1];  %function for derivative of PD matrices
        end
        ss_dist_min = [0.1 10];
    elseif gs_num == 1
        Theta2 = 0;
        d_Theta = {[-17.46 17.46], 0}; %Droumax and Droumin  % [-17.46 17.46]
        if XY_quadterm == 0 % without second order term
            FthetaNum = [1 1];
            F_theta = @(x) [1 x(1)]; %function for PD matrices
            d_F_theta = @(x) [0 1];  %function for derivative of PD matrices            
        else % with second order term
            FthetaNum = [1 2];
            F_theta = @(x) [1 x(1) x(1)^2]; %function for PD matrices
            d_F_theta = @(x) [0 1 2*x(1)];
        end   
        ss_dist_min = 0.1;
    end     
   % [1.464 13.642]
end

%% Hinf controller design
for active = 1:0
    theta_eval = [6.8213;3000]; % maDot:0.5, N:3000
    ny = size(Gasym.C2,1);
    nu = size(Gasym.B2,2);
    Ga = AugPltEval(Gasym,theta_eval);
    a = Ga.A;
    b = [Ga.B1 Ga.B2];
    c = [Ga.C1;Ga.C2];
    d = [Ga.D11 Ga.D12;
        Ga.D21 Ga.D22];
    ga = ss(a,b,c,d);    
    [Khinf,CL,gam,info] = hinfsyn(ga,ny,nu,'method','lmi');
%             tf(CL)
    assignin('base','Khinf',Khinf);
    assignin('base','Gam',gam);
gam
end

% prepare the data for output
theta_info.Theta1 = Theta1;
theta_info.Theta2 = Theta2;
theta_info.d_Theta = d_Theta;
theta_info.gs_num = gs_num;
theta_info.ss_dist_min = ss_dist_min;

theta_fcn_info.F_theta = F_theta;
theta_fcn_info.d_F_theta = d_F_theta;
theta_fcn_info.FthetaNum = FthetaNum;

if ~exist('num_gain')
    num_gain = 1;
end
