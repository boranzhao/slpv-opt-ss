function Gam = gamma_cal_slpv(tau, Gasym, Ga_ctmat, theta_info, design_para,opt_lmilab)
regnum1 = size(theta_info.Theta1,1);
regnum2 = size(theta_info.Theta2,1);
tau1 = tau(1:(regnum1-1)*2);
tau2 = tau((regnum1-1)*2+1:end);

if regnum1 == 2
    theta_info.Theta1 = [theta_info.Theta1(1,1) tau1(2);tau1(1) theta_info.Theta1(2,2)];
elseif regnum1 == 3
    theta_info.Theta1 = [theta_info.Theta1(1,1) tau1(2);tau1(1) tau1(4);tau1(3) theta_info.Theta1(3,2)];
end
if regnum2 == 2
    theta_info.Theta2 = [theta_info.Theta2(1,1) tau2(2);tau2(1) theta_info.Theta2(2,2)];
elseif regnum2 == 3
    theta_info.Theta2 = [theta_info.Theta2(1,1) tau2(2);tau2(1) tau2(4);tau2(3) theta_info.Theta2(3,2)];
end 
if design_para.LPVdesign_method == 1
    [OptRst,sol] = SLPV_Basic_Affine_MC(Gasym,Ga_ctmat,theta_info,design_para,opt_lmilab);
elseif  design_para.LPVdesign_method == 2                    
    [OptRst,sol] = SLPV_Proj_Affine_MC(Gasym,Ga_ctmat,theta_info,design_para,opt_lmilab);
end
Gam = OptRst.Gam; %% changed 