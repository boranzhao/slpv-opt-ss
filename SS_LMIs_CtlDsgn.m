function [lmisys,lminum] = SS_LMIs_CtlDsgn(lmisys0,X,Y,Ftheta,jr,design_para,lminum,X_k,Y_k)
%% Code for generating LMIs based on LMI Lab at switching surfaces given optimization variables 
XY_PD = design_para.XY_PD;
eps_eq = design_para.eps_eq;
ss_cond = design_para.ss_cond;
mu = 1e6*eps_eq;
if nargin < 8
    X_k = [];
    Y_k = [];
end

j = jr(1); r = jr(2);
setlmis(lmisys0)
switch XY_PD 
    case 0 %for constant LF, Lyapunov functions have to be same for all subsets. 
    case 1  % PD X, constant Y  
         % X(j) >= X(r)
        lminum = lminum + 1;
        for id_Ftheta = 1:length(Ftheta)
            lmiterm([lminum 1 1 X(id_Ftheta,j)],-Ftheta(id_Ftheta), 1);
            lmiterm([lminum 1 1 X(id_Ftheta,r)],Ftheta(id_Ftheta), 1);
        end 
    case 2 % PD Y, so constant X should be equal at SS
        % Y(j) <= Y(r)
        lminum = lminum + 1;
        for id_Ftheta = 1:length(Ftheta)
            lmiterm([lminum 1 1 Y(id_Ftheta,j)],Ftheta(id_Ftheta), 1);
            lmiterm([lminum 1 1 Y(id_Ftheta,r)],-Ftheta(id_Ftheta), 1);
        end      
    case 3
%          Consts = [Consts [eps_eq*In Y-Y1; Y-Y1 eps_eq*In]>=0 X>=X1];
            if ss_cond== 1
                if design_para.ss_improve == 0 || isempty(Y_k)
    %                 % Yj == Yr, 
    %                 lminum = lminum + 1;
    %                 lmiterm([lminum 1 1 0],-eps_eq);
    %                 lmiterm([lminum 2 2 0],-eps_eq);
    %                 for id_Ftheta = 1:length(Ftheta)
    %                     lmiterm([lminum 2 1 Y(id_Ftheta,j)],-Ftheta(id_Ftheta), 1);
    %                     lmiterm([lminum 2 1 Y(id_Ftheta,r)],Ftheta(id_Ftheta), 1);                                   
    %                 end    
                    % Yj <= Yr,  Yr-Yj <= eps*I
                    if design_para.ss_cvx.method == 2
                        lminum = lminum + 1;
                        for id_Ftheta = 1:length(Ftheta)
                            lmiterm([lminum 1 1 Y(id_Ftheta,j)],Ftheta(id_Ftheta),1);
                            lmiterm([lminum 1 1 Y(id_Ftheta,r)],-Ftheta(id_Ftheta),1);
                        end

                        lminum = lminum + 1;
                        for id_Ftheta = 1:length(Ftheta)
                            lmiterm([lminum 1 1 Y(id_Ftheta,r)],Ftheta(id_Ftheta),1);
                            lmiterm([lminum 1 1 Y(id_Ftheta,j)],-Ftheta(id_Ftheta),1);
                        end
                        lmiterm([lminum 1 1 0],-eps_eq);
                    elseif design_para.ss_cvx.method == 3 % relax original nonconvex constraint by adding an extra constraint inv(Y)>=eps*I 
                        lminum = lminum + 1;
                        for id_Ftheta = 1:length(Ftheta)
                            lmiterm([lminum 1 1 Y(id_Ftheta,j)],Ftheta(id_Ftheta),1);
                            lmiterm([lminum 1 1 Y(id_Ftheta,r)],-Ftheta(id_Ftheta),1);
                        end
                        
                        lminum = lminum + 1;
                        for id_Ftheta = 1:length(Ftheta)
                            lmiterm([lminum 1 1 Y(id_Ftheta,r)],Ftheta(id_Ftheta),1);
                        end
                        lmiterm([lminum 1 1 0],-1/design_para.ss_cvx.eps);
                        
                        lminum = lminum + 1;
                        for id_Ftheta = 1:length(Ftheta)
                            lmiterm([-lminum 1 1 X(id_Ftheta,j)],Ftheta(id_Ftheta), 1);
                            lmiterm([-lminum 1 1 X(id_Ftheta,r)],-Ftheta(id_Ftheta), 1);
                            lmiterm([-lminum 2 2 Y(id_Ftheta,j)],Ftheta(id_Ftheta), 1);
                        end   
                        lmiterm([-lminum 1 1 0],design_para.ss_cvx.eps); % to be safe  
                        lmiterm([-lminum 2 1 0],1); % to be safe                    
                    end
                    
                    if design_para.ss_cvx.method ~= 3 
                        % Xj >= Xr+mu*I;
                        lminum = lminum + 1;  
                        for id_Ftheta = 1:length(Ftheta)
                            lmiterm([lminum 1 1 X(id_Ftheta,j)],-Ftheta(id_Ftheta), 1);
                            lmiterm([lminum 1 1 X(id_Ftheta,r)],Ftheta(id_Ftheta), 1);
                        end   
    %                     lmiterm([lminum 1 1 0],eps_eq); % to be safe
                        if design_para.ss_cvx.method == 2
                            lmiterm([lminum 1 1 0],eps_eq);
                        end                    
                        % Yr >= I to guarantee Xj - Yj^-1 >= Xr - Yr^-1
    %                     lminum = lminum + 1;  
    %                     for id_Ftheta = 1:length(Ftheta)
    %                         lmiterm([lminum 1 1 Y(id_Ftheta,r)],-Ftheta(id_Ftheta), 1);
    %                     end    
    %                     lmiterm([lminum 1 1 0],sqrt(eps_eq/mu)); 
                    end
                elseif design_para.ss_improve == 1 && ~isempty(Y_k)
                    % Yj <= Yr, 
                    lminum = lminum + 1;
                    for id_Ftheta = 1:length(Ftheta)
                        lmiterm([lminum 1 1 Y(id_Ftheta,j)],Ftheta(id_Ftheta),1);
                        lmiterm([lminum 1 1 Y(id_Ftheta,r)],-Ftheta(id_Ftheta),1);                
                    end
                    lmiterm([lminum 1 1 0],-eps_eq); % for relaxing
                    % Xj-Yj^-1>= Xr-Yr^-1;
                    n = size(Y_k,1);
                    Y_r_k = zeros(n);
                    for id_Ftheta = 1:length(Ftheta)
                        Y_r_k = Y_r_k + Ftheta(id_Ftheta)*Y_k(:,:,id_Ftheta,r);
                    end    
                    lminum = lminum + 1;                  
                    for id_Ftheta = 1:length(Ftheta)
                        lmiterm([-lminum 1 1 X(id_Ftheta,j)],Ftheta(id_Ftheta)*Y_r_k,Y_r_k); 
                        lmiterm([-lminum 1 1 X(id_Ftheta,r)],-Ftheta(id_Ftheta)*Y_r_k,Y_r_k);             
                        lmiterm([-lminum 1 1 Y(id_Ftheta,r)],-Ftheta(id_Ftheta),1);                  
                        lmiterm([-lminum 2 2 Y(id_Ftheta,j)],Ftheta(id_Ftheta),1); 
                    end            
                    lmiterm([-lminum 1 1 0],2*Y_r_k);             
                    lmiterm([-lminum 2 1 0],Y_r_k);   
                end                    
            elseif ss_cond == 2                                
%          Consts = [Consts [eps_eq*In X-X1; X-X1 eps_eq*In]>=0 Y<=Y1];
%                 % Xj == Xr,                 
%                 lminum = lminum + 1;
%                 lmiterm([lminum 1 1 0],-eps_eq);
%                 lmiterm([lminum 2 2 0],-eps_eq);
%                 for id_Ftheta = 1:length(Ftheta)
%                     lmiterm([lminum 2 1 X(id_Ftheta,j)],-Ftheta(id_Ftheta), 1);
%                     lmiterm([lminum 2 1 X(id_Ftheta,r)],Ftheta(id_Ftheta), 1);                                  
%                 end  
                if design_para.ss_improve == 0 || isempty(X_k)
                    if design_para.ss_cvx.method == 2
                        % Xj >= Xr,  Xj-Xr <= eps*I
                        lminum = lminum + 1;
                        for id_Ftheta = 1:length(Ftheta)
                            lmiterm([lminum 1 1 X(id_Ftheta,r)],Ftheta(id_Ftheta),1);
                            lmiterm([lminum 1 1 X(id_Ftheta,j)],-Ftheta(id_Ftheta),1);
                        end

                        lminum = lminum + 1;
                        for id_Ftheta = 1:length(Ftheta)
                            lmiterm([lminum 1 1 X(id_Ftheta,j)],Ftheta(id_Ftheta),1);
                            lmiterm([lminum 1 1 X(id_Ftheta,r)],-Ftheta(id_Ftheta),1);
                        end
                        lmiterm([lminum 1 1 0],-eps_eq);
                    elseif design_para.ss_cvx.method == 3
                        lminum = lminum + 1;
                        for id_Ftheta = 1:length(Ftheta)
                            lmiterm([lminum 1 1 X(id_Ftheta,r)],Ftheta(id_Ftheta),1);
                            lmiterm([lminum 1 1 X(id_Ftheta,j)],-Ftheta(id_Ftheta),1);
                        end
                        
                        lminum = lminum + 1;
                        for id_Ftheta = 1:length(Ftheta)
                            lmiterm([lminum 1 1 X(id_Ftheta,j)],Ftheta(id_Ftheta),1);
                        end
                        lmiterm([lminum 1 1 0],-1/design_para.ss_cvx.eps);
                        
                        lminum = lminum + 1;
                        for id_Ftheta = 1:length(Ftheta)
                            lmiterm([-lminum 1 1 Y(id_Ftheta,r)],Ftheta(id_Ftheta), 1);
                            lmiterm([-lminum 1 1 Y(id_Ftheta,j)],-Ftheta(id_Ftheta), 1);
                            lmiterm([-lminum 2 2 X(id_Ftheta,r)],Ftheta(id_Ftheta), 1);
                        end   
                        lmiterm([-lminum 1 1 0],design_para.ss_cvx.eps); % to be safe  
                        lmiterm([-lminum 2 1 0],1); % to be safe    
                    end
                    if design_para.ss_cvx.method ~= 3
                        % Yj < = Yr -eps_eq*I
                        lminum = lminum + 1;
                        for id_Ftheta = 1:length(Ftheta)
                            lmiterm([lminum 1 1 Y(id_Ftheta,j)],Ftheta(id_Ftheta), 1);
                            lmiterm([lminum 1 1 Y(id_Ftheta,r)],-Ftheta(id_Ftheta), 1);
                        end 
                        if design_para.ss_cvx.method == 2
                            lmiterm([lminum 1 1 0],eps_eq);
                        end
                    end
                elseif design_para.ss_improve == 1 && ~isempty(X_k)
                    % Xj >= Xr, 
                    lminum = lminum + 1;
                    for id_Ftheta = 1:length(Ftheta)
                        lmiterm([lminum 1 1 X(id_Ftheta,r)],Ftheta(id_Ftheta),1);
                        lmiterm([lminum 1 1 X(id_Ftheta,j)],-Ftheta(id_Ftheta),1);
                    end
                    lmiterm([lminum 1 1 0],-eps_eq); % for relaxing

%                     % Yj- Xj^-1 <= Yr - Xr^-1;
                    n = size(X_k,1);
                    X_j_k = zeros(n);
                    for id_Ftheta = 1:length(Ftheta)
                        X_j_k = X_j_k + Ftheta(id_Ftheta)*X_k(:,:,id_Ftheta,j);
                    end    
                    lminum = lminum + 1;                  
                    for id_Ftheta = 1:length(Ftheta)
                        lmiterm([-lminum 1 1 Y(id_Ftheta,r)],Ftheta(id_Ftheta)*X_j_k,X_j_k); 
                        lmiterm([-lminum 1 1 Y(id_Ftheta,j)],-Ftheta(id_Ftheta)*X_j_k,X_j_k);             
                        lmiterm([-lminum 1 1 X(id_Ftheta,j)],-Ftheta(id_Ftheta),1);                  
                        lmiterm([-lminum 2 2 X(id_Ftheta,r)],Ftheta(id_Ftheta),1); 
                    end            
                    lmiterm([-lminum 1 1 0],2*X_j_k);             
                    lmiterm([-lminum 2 1 0],X_j_k);   
%                    just for test
%                     lmiterm([-lminum 1 1 0],eps_eq*10);  
%                     lmiterm([-lminum 2 2 0],eps_eq*10); 
                end                          
            end 
end % switch
lmisys = getlmis;
end


