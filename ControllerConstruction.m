%% controller construction

% Note that D11,D12 and D21 have to be parameter independent; otherwise,
% the following process has to be implemented during each sampling interval
% compute Dk (nu*ny) s.t. sig(D11+D12DkD21) < gam for each subset
gam = OptRst.gam;
Dk_lmivar = zeros(1,regnum);
Dk = zeros(nu,ny,regnum);
opt_lmilab(1)= 1e-5;    % relative accuary on the optimal value
opt_lmilab(2)= 600;     % Number of Iteration
opt_lmilab(4) = 10;     % J, the code terminates when the objective has not decreased by more than the desired relative accuracy during the last J iterations
opt_lmilab(5)=  1;
for regid = 1:regnum
    setlmis([]);
    Dk_lmivar(regid) = lmivar(2,[nu ny]);
    lminum = 1;
    
    lmiterm([-lminum 1 1 0],gam(regid));
    lmiterm([-lminum 2 1 Dk_lmivar(regid)],Gasym.D12, Gasym.D21);
    lmiterm([-lminum 2 1 0],Gasym.D11);
    lmiterm([-lminum 2 2 0],gam(regid));   
    
        
    lmisys = getlmis;
    [~,Dk(:,:,regid)] = feasp(lmisys,opt_lmilab);    
end




