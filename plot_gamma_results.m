%% for numerical example
% LPV: 17.96
% for max gamma_j
% ss                gamma
% this [-0.1 0.1]        [17.06, 17.06]
% [0.6 0.8]         [16.17  16.17]
% this [0.5 0.7]    [16.35  16.35]
% [0.5 0.8]         [16.28  16.28]
% optimal ss:  [0.80096 0.81214]   [13.72 13.72]  

% for weighted average of 
% ss                gamma
% [-0.1 0.1]        [14.625, 17.95]
% [-0.1 0.1]        [16.4, 17.23] % [3.87 -1]% improvement over the traditional one=
% this: [-0.1 0.1]        [15.83, 17.40] % [7.2 -2]% improvement than the traditional one
% [0.6 0.8]]        [14.75  17.96]
% [0.6 0.8]         [15.56 16.49] % [3.8 -2]% improvement than the traditional one  
% [0.6 0.8]         [15.56 16.49] % [4.8 -3]% worse than the traditonal one
% this [0.5 0.7]         [15.65 16.68] % [4.3 -2]% worse than th

% [0.5 0.8]         [15.68 16.61] % [3.7 -2]% better than th
% optimal SS: [0.80 0.81]       [13.14 14.92]  % use LPV as the bound
% optimal SS: [0.8095 0.8197]   [13.36 13.9944] % [2.6 -2]% improvement
% compared to optimal traditonal SLPV

close all; 
linewd = 1;
Gam_NS = 17.959;
FigIndex = 3;

if FigIndex == 1
    SLPV0.gam = [17.059 17.059];
    SLPV1.gam = [15.83, 17.40] ;
    SLPV0.Theta1 = [-1 0.1;-0.1 1];
    SLPV1.Theta1 = [-1 0.1;-0.1 1];
elseif FigIndex == 2
    SLPV0.gam = [16.35 16.35];
    SLPV1.gam = [15.65, 16.68] ;
    SLPV0.Theta1 = [-1 0.8;0.5 1];   
    SLPV1.Theta1 = [-1 0.8;0.5 1];   
elseif FigIndex == 3
    SLPV0.gam = [13.72, 13.72];
    SLPV1.gam = [13.14, 14.92];    
    SLPV0.Theta1 = [-1 0.81214;0.80096 1];
    SLPV1.Theta1 = [-1 0.81;0.80 1];  
end
% OptSLPV_optSS.gam = [13.137, 14.920];
% theta_info_optSS.Theta1 = [-1 0.810;0.800 1];

height = 0.2;
legend_str = cell(1,1);
h =[];
id = 1;
regnum = length(SLPV1.gam);
h(id) = plot([SLPV0.Theta1(1,1) SLPV0.Theta1(end,end)],[Gam_NS Gam_NS],'k-.','Linewidth',linewd*2);
legend_str{id} = 'LPV';
id = id + 1;
hold on;
 
for i=1:regnum
    h(id) = plot([SLPV0.Theta1(i,1) SLPV0.Theta1(i,2)],[SLPV0.gam(i) SLPV0.gam(i)],'m:','Linewidth',linewd*2);
end
plot([SLPV0.Theta1(1,2) SLPV0.Theta1(1,2)],[SLPV0.gam(1)-height SLPV0.gam(1)+height],'m-.','Linewidth',linewd)
plot([SLPV0.Theta1(2,1) SLPV0.Theta1(2,1)] ,[SLPV0.gam(2)-height SLPV0.gam(2)+height],'m-.','Linewidth',linewd), 
if FigIndex == 1
    legend_str{id} = 'SLPV 0(a)' ;
elseif FigIndex == 2
    legend_str{id} = 'SLPV 0(b)' ;
else
    legend_str{id} = 'SLPV 0(c)' ;
end
id = id + 1;
 
% switching LPV controller with new cost function (max gamma_j)
% and optimized switching surfaces
 

for i=1:regnum
    h(id) = plot([SLPV0.Theta1(i,1) SLPV0.Theta1(i,2)],[SLPV1.gam(i) SLPV1.gam(i)],'r-','Linewidth',linewd*1);
end
plot([SLPV1.Theta1(1,2) SLPV1.Theta1(1,2)],[SLPV1.gam(1)-height SLPV1.gam(1)+height],'r-.','Linewidth',linewd)
plot([SLPV1.Theta1(2,1) SLPV1.Theta1(2,1)] ,[SLPV1.gam(2)-height SLPV1.gam(2)+height],'r-.','Linewidth',linewd),
if FigIndex == 1
    legend_str{id} = 'SLPV 1(a)' ;
elseif FigIndex == 2
    legend_str{id} = 'SLPV 1(b)' ;
else
    legend_str{id} = 'SLPV 1(c)' ;
end
id = id + 1;


% another proposed SLPV controller with only 2perc worse than the traditonal controller
SLPV1.gam = [13.36, 13.99];    
SLPV1.Theta1 = [-1 0.8203;0.8101 1];  
for i=1:regnum
    h(id) = plot([SLPV0.Theta1(i,1) SLPV0.Theta1(i,2)],[SLPV1.gam(i) SLPV1.gam(i)],'b--','Linewidth',linewd*1.5);
end
plot([SLPV1.Theta1(1,2) SLPV1.Theta1(1,2)],[SLPV1.gam(1)-height SLPV1.gam(1)+height],'b-','Linewidth',linewd)
plot([SLPV1.Theta1(2,1) SLPV1.Theta1(2,1)] ,[SLPV1.gam(2)-height SLPV1.gam(2)+height],'b-.','Linewidth',linewd),
legend_str{id} = 'SLPV 1(d)' ;

% optimized switching surfaces and iteration on ss condition 
% for i=1:regnum
%     h4 = plot([theta_info_optSS.Theta1(i,1) theta_info_optSS.Theta1(i,2)],[OptSLPV_optSS.gam(i) OptSLPV_optSS.gam(i)],'r','Linewidth',linewd*2);
% end
% plot([theta_info_optSS.Theta1(1,2) theta_info_optSS.Theta1(1,2)],[OptSLPV_optSS.gam(1)-height OptSLPV_optSS.gam(1)+height],'r-.','Linewidth',linewd)
% plot([theta_info_optSS.Theta1(2,1) theta_info_optSS.Theta1(2,1)] ,[OptSLPV_optSS.gam(2)-height OptSLPV_optSS.gam(2)+height],'r-.','Linewidth',linewd), 
%  
% legend_str{id} = 'SLPV 1(c)'; % 'SLPV with ALLC and optimized SSs' ;
% id = id+1;
% switching LPV controller with new cost funct
 
legend(h,legend_str,'Location','best');
xlabel('\theta')
ylabel('\gamma-value');
% ylim([15 18]);
goodplot([7 4],[0.01 0.2]);
%%
pause;
%%
if FigIndex == 1  
    print -painters -dpdf -r150 gamCompPlt6_1.pdf
elseif FigIndex == 2
    print -painters -dpdf -r150 gamCompPlt6_2.pdf
elseif FigIndex == 3
    print -painters -dpdf -r150 gamCompPlt6_3.pdf
end



