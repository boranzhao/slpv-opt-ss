%% Plot data
% clear;
close all;
clear;
clear legend
%clc;
DataSel = 1;
Hinf = 'Hinf_3000_jnl_nojump.mat';
LPV = 'LPV_XY3_Proj_jnl1.mat';
SLPV_heuSS = 'SLPV_heuSS_3subs_XY3_Proj_jnl_nojump.mat';
SLPV_heuSS_max = 'SLPV_heuSS_3subs_XY3_Proj_max.mat';
SLPV_optSS = 'SLPV_optSS_3subs_XY3_Proj_jnl_nojump.mat';

if DataSel == 1 % sinusoidal
    if exist(Hinf,'file') == 2
        load(Hinf,'Data');
        DataHinf = Data;
    end
    if exist(LPV,'file') == 2
        load(LPV,'Data','OptRst');
        DataLPV = Data;
        OptLPV = OptRst;
    end    
    if exist(SLPV_heuSS,'file') == 2
        load(SLPV_heuSS,'Data','theta_info','OptRst');
        DataSLPV_heuSS = Data;
        OptSLPV_heuSS = OptRst;
        theta_info_heuSS = theta_info;
    end
    if exist(SLPV_heuSS_max,'file') == 2
        load(SLPV_heuSS_max,'Data','theta_info','OptRst');
        DataSLPV_heuSS_max = Data;
        OptSLPV_heuSS_max = OptRst;
        theta_info_heuSS_max = theta_info;
    end
    if exist(SLPV_optSS,'file') == 2
        load(SLPV_optSS,'Data','theta_info','OptRst');
        DataSLPV_optSS = Data;
        OptSLPV_optSS = OptRst;
        theta_info_optSS = theta_info;
    end   
elseif DataSel == 2 % 
    load('LPV_PDY_Step.mat','LPVData' );
    load('SLPV_PDY_Step_IgnoreUncert.mat','SLPVData');
    SLPVData_IgUn = SLPVData;
    loclearad('SLPV_PDY_Step.mat','SLPVData');
    load('LPV_PDY_Step.mat','w' );
end
%% 
% line style:
% Hinf: dashed line, thin, black
% LPV: dash-dot line, thick, magenta
% SLPV_heuSS: solid line, thin, blue
% SLPV_optSS: solid line, thick, red
figure;
subplot(2,1,1)
plot(DataLPV.r.Time,DataLPV.maDot.Data*100,'k','LineWidth',1.5);
xlabel('Time (s)')
ylabel('Air flow (%)');
axis([0 60 0 110]);
goodplot
subplot(2,1,2)
plot(DataLPV.r.Time,DataLPV.N.Data,'k','LineWidth',1.5);
xlabel('Time (s)')
ylabel('Engine Speed (RPM)')
axis([0 60 0 6200]);
goodplot([7 5]);
% print -painters -dpdf -r150 maN_profile_jnl.pdf

%% XY plot
figure;
% subplot(1,2,1)
% plot(DataLPV.N.Data,DataLPV.maDot.Data*100,'k','LineWidth',1.5);
% xlabel('Engine Speed (RPM)');
% ylabel('Air flow (%)');
% axis([0 6000 0 100]);
% goodplot
% figure;
% % subplot(1,2,2);
% plot([800 6000],theta_info_optSS.Theta1(2,1)*[1, 1], 'r--','LineWidth',1);
% hold on;
% plot([800 6000],theta_info_optSS.Theta1(1,2)*[1, 1],'r--','LineWidth',1);
% plot([800 6000],theta_info_optSS.Theta1(3,1)*[1, 1], 'r--','LineWidth',1);
% plot([800 6000],theta_info_optSS.Theta1(2,1)*[2, 2],'r--','LineWidth',1);
% hold on;
% plot(DataLPV.N.Data,DataLPV.Tinv.Data,'k','LineWidth',1.5);
% xlabel('Engine Speed (RPM)');
% ylabel('Inverse of time delay (s^{-1})');
% axis([0 6000 0 15]);
% hold on;
gain = 1000;
plot(DataLPV.maDot.Data*100,DataLPV.N.Data/gain,'k','LineWidth',1.5);
ylabel('Engine Speed (\times1000 RPM)');
xlabel('Air flow (%)');
axis([ 0 100 500/gain 6000/gain]);
goodplot([4 4]);
% print -painters -dpdf -r150 maN_XYplot_jnl.pdf
figure;
% subplot(1,2,2);
plot(theta_info_optSS.Theta1(2,1)*[1, 1],[800 6000]/gain,'r-.','LineWidth',1);
hold on;
plot(theta_info_optSS.Theta1(1,2)*[1, 1],[800 6000]/gain,'r-.','LineWidth',1);
plot(theta_info_optSS.Theta1(3,1)*[1, 1],[800 6000]/gain, 'r-.','LineWidth',1);
plot(theta_info_optSS.Theta1(2,2)*[1, 1],[800 6000]/gain,'r-.','LineWidth',1);
plot(theta_info_optSS.Theta1(1,1)*[1, 1],[800 6000]/gain, 'b--', 'LineWidth',1);
plot(theta_info_optSS.Theta1(3,2)*[1, 1],[800 6000]/gain, 'b--', 'LineWidth',1);
plot(DataLPV.Tinv.Data,DataLPV.N.Data/gain,'k','LineWidth',1.5);
ylabel('\theta_2 (\times1000 RPM)');
xlabel('\theta_1 (s^{-1})');
axis([1 14 500/gain 6000/gain ]);
hold on;
goodplot([4 4]);
% print -painters -dpdf -r150 TinvN_XYplot_jnl.pdf
%%
t = 0:0.001:60;
legend_txt = cell(1,1); i=1;
legend_txt2 = cell(1,1);j=1;
h_u = figure;
hold on;
h = figure; 
% h_sig = figure;
subplot(2,1,1);

plot(DataLPV.r.Time,DataLPV.r.Data,'k--','Linewidth',1); % the reference line
legend_txt{i}='Reference'; i = i+1;
hold on;
subplot(2,1,2);
hold on;
if exist ('DataHinf','var')
   subplot(2,1,1)
   plot(DataHinf.r.Time,DataHinf.y.Data,':','color',[0 0.5 1],'Linewidth',1.5);
   legend_txt{i} = 'H_\infty'; 
   i=i+1;
   err = (DataHinf.y.Data-DataHinf.r.Data);
   err = interp1(DataHinf.r.Time,err,t);
%    figure;
%    plot(DataHinf.r.Time,err,t,err_s);
   Hinf_err = sqrt(sum(err.^2)/length(err))
   figure(h_u);
   plot(DataHinf.r.Time,DataHinf.u.Data,':','color',[0 0.5 1],'Linewidth',1.5);
end
if exist('DataLPV','var')
    figure(h);
    subplot(2,1,1);
    plot(DataLPV.r.Time,DataLPV.y.Data,'m','Linewidth',1);  
    legend_txt{i} = 'LPV'; 
    i=i+1;
    err = (DataLPV.y.Data-DataLPV.r.Data);
    err = interp1(DataLPV.r.Time,err,t,'linear');
%     figure;
%     plot(DataLPV.r.Time,err);hold on;
%     plot(t,err_s);
    LPV_err = sqrt(sum(err.^2)/length(err))
    figure(h_u);
    plot(DataLPV.r.Time,DataLPV.u.Data,'m','Linewidth',1);  
end
if exist('DataSLPV_heuSS','var')
    figure(h);
    subplot(2,1,1)
    plot(DataSLPV_heuSS.r.Time,DataSLPV_heuSS.y.Data,'b-.','Linewidth',1.5);  
    legend_txt{i} = 'SLPV 1'; % 'SLPV with ALLCs';
    i=i+1;
    subplot(2,1,2)
    plot(DataSLPV_heuSS.r.Time,DataSLPV_heuSS.CtrlSel.Data,'b-.','Linewidth',1.5);  
    legend_txt2{j} = 'SLPV 1';
    j = j+1;
    err = (DataSLPV_heuSS.y.Data-DataSLPV_heuSS.r.Data);
    err = interp1(DataSLPV_heuSS.r.Time,err,t,'linear');
%     figure;
%     plot(DataSLPV_heuSS.r.Time,err);hold on;
%     plot(t,err_s);
    SLPV_heu_err = sqrt(sum(err.^2)/length(err))   
    figure(h_u);
    plot(DataSLPV_heuSS.r.Time,DataSLPV_heuSS.u.Data,'b-.','Linewidth',1.5);  
end
if exist('DataSLPV_optSS','var')
    figure(h);
    subplot(2,1,1)
    plot(DataSLPV_optSS.r.Time,DataSLPV_optSS.y.Data,'r-','Linewidth',1);
    legend_txt{i} = 'SLPV 2'; %'SLPV with ALLCs and optimized SSs';
    i=i+1;
    subplot(2,1,2)
    plot(DataSLPV_optSS.r.Time,DataSLPV_optSS.CtrlSel.Data,'r-','Linewidth',1);  
    legend_txt2{j} = 'SLPV 2';
    j = j+1;
    err = (DataSLPV_optSS.y.Data-DataSLPV_optSS.r.Data);
    err = interp1(DataSLPV_optSS.r.Time,err,t,'linear');
%     figure;
%     plot(DataSLPV_optSS.r.Time,err);hold on;
%     plot(t,err_s);
    SLPV_opt_err = sqrt(sum(err.^2)/length(err))
    figure(h_u);
    plot(DataSLPV_optSS.r.Time,DataSLPV_optSS.u.Data,'r','Linewidth',1);
end

figure(h);
subplot(2,1,1)
% plot([0 60],[1.1 1.1],'k:'); % the reference line
xlabel('Time (s)');
ylabel('Equivalence ratio');
legend(legend_txt(2:end),'FontSize',10,'Location','best','Orientation','horizontal'); %
legend('boxoff')
ylim([0.7 1.3]);
% ylim([-1.2 1.2]);
xlim([20 60]);
grid on;
% rectangle('Position', [20 0.85 40 0.3],'EdgeColor','k','LineStyle','-.','Linewidth',1);
goodplot([7,6])


subplot(2,1,2)
xlabel('Time (s)');
ylabel('Switching signal');
legend(legend_txt2,'FontSize',10,'Location','best'); %,'Orientation','horizontal'
legend('boxoff')
goodplot([7,6])

% print -painters -dpdf -r150 dist_rej_jnl.pdf

figure(h_u)
xlabel('Time (s)');
ylabel('Control input');
legend(legend_txt(2:end),'FontSize',10,'Location','best');
pause;
%% compare the subset gamma value yielded by different controllers
height = 0.08;
figure;
linewd = 1;
id = 1;
Gam_NS = OptLPV.Gam;
regnum = length(OptSLPV_heuSS.gam);
h(id) = plot([theta_info.Theta1(1,1) theta_info.Theta1(end,end)],[Gam_NS Gam_NS],'k-.','Linewidth',linewd*1.5);
legend_str{id} = 'LPV';
hold on;

id=id+1;
legend_str{id} = 'SLPV 0'; %  'using the traditional approach' ;
for i=1:regnum
    h(id) = plot([theta_info_heuSS.Theta1(i,1) theta_info_heuSS.Theta1(i,2)],[OptSLPV_heuSS_max.gam(i) OptSLPV_heuSS_max.gam(i)],'m:','Linewidth',linewd*1.5);
end
plot([theta_info_heuSS.Theta1(1,2) theta_info_heuSS.Theta1(1,2)],[OptSLPV_heuSS_max.gam(1)-height OptSLPV_heuSS_max.gam(1)+height],'m-.','Linewidth',linewd)
plot([theta_info_heuSS.Theta1(2,1) theta_info_heuSS.Theta1(2,1)] ,[OptSLPV_heuSS_max.gam(2)-height OptSLPV_heuSS_max.gam(2)+height],'m-.','Linewidth',linewd),
plot([theta_info_heuSS.Theta1(2,2) theta_info_heuSS.Theta1(2,2)],[OptSLPV_heuSS_max.gam(2)-height OptSLPV_heuSS_max.gam(2)+height],'m-.','Linewidth',linewd)
plot([theta_info_heuSS.Theta1(3,1) theta_info_heuSS.Theta1(3,1)] ,[OptSLPV_heuSS_max.gam(3)-height OptSLPV_heuSS_max.gam(3)+height],'m-.','Linewidth',linewd), 

% SLPV with improved local perforemance
id=id+1;
legend_str{id} = 'SLPV 1'; 
for i=1:regnum
    h(id) = plot([theta_info_heuSS.Theta1(i,1) theta_info_heuSS.Theta1(i,2)],[OptSLPV_heuSS.gam(i) OptSLPV_heuSS.gam(i)],'b--','Linewidth',linewd*1.5);
end
plot([theta_info_heuSS.Theta1(1,2) theta_info_heuSS.Theta1(1,2)],[OptSLPV_heuSS.gam(1)-height OptSLPV_heuSS.gam(1)+height],'b-.','Linewidth',linewd)
plot([theta_info_heuSS.Theta1(2,1) theta_info_heuSS.Theta1(2,1)] ,[OptSLPV_heuSS.gam(2)-height OptSLPV_heuSS.gam(2)+height],'b-.','Linewidth',linewd), 
plot([theta_info_heuSS.Theta1(2,2) theta_info_heuSS.Theta1(2,2)],[OptSLPV_heuSS.gam(2)-height OptSLPV_heuSS.gam(2)+height],'b-.','Linewidth',linewd)
plot([theta_info_heuSS.Theta1(3,1) theta_info_heuSS.Theta1(3,1)] ,[OptSLPV_heuSS.gam(3)-height OptSLPV_heuSS.gam(3)+height],'b-.','Linewidth',linewd), 

% SLPV with improved local performance and optimized SSs
id=id+1;
legend_str{id} =  'SLPV 2'; %'SLPV with ALLC and optimized SSs' ;
for i=1:regnum
    h(id) = plot([theta_info_optSS.Theta1(i,1) theta_info_optSS.Theta1(i,2)],[OptSLPV_optSS.gam(i) OptSLPV_optSS.gam(i)],'r','Linewidth',linewd*2);
end
plot([theta_info_optSS.Theta1(1,2) theta_info_optSS.Theta1(1,2)],[OptSLPV_optSS.gam(1)-height OptSLPV_optSS.gam(1)+height],'r-.','Linewidth',linewd)
plot([theta_info_optSS.Theta1(2,1) theta_info_optSS.Theta1(2,1)] ,[OptSLPV_optSS.gam(2)-height OptSLPV_optSS.gam(2)+height],'r-.','Linewidth',linewd), 
plot([theta_info_optSS.Theta1(2,2) theta_info_optSS.Theta1(2,2)],[OptSLPV_optSS.gam(2)-height OptSLPV_optSS.gam(2)+height],'r-.','Linewidth',linewd)
plot([theta_info_optSS.Theta1(3,1) theta_info_optSS.Theta1(3,1)] ,[OptSLPV_optSS.gam(3)-height OptSLPV_optSS.gam(3)+height],'r-.','Linewidth',linewd), 

% switching LPV controller with new cost funct
% optimized switching surfaces and iteration on ss condition

legend(h, legend_str,'Location','east');%'Orientation','horizontal'
% legend('boxoff')
xlabel('Scheduling parameter \theta_1')
ylabel('\gamma value');
goodplot;
%%
% print -painters -dpdf -r150 gam_comp_jnl.pdf





