% Emily Duan
% DISPLAY RESULTS FROM OPTIMIZATION
clc
clear all
close all

% LOAD RESULTS FROM OPTIMIZATION AND EXTRACT FAM BUNDLE DESIGN INFO
load('maxFbkh_maxdeltalm_5000gen_rect_lateral_mar19_6x12x1_0.3E_0.5P_0.3M_0.1667h0_trial1_C01_260.5543e3.mat')
POPULATION_TRIAL1 = POPULATION(num_gen+1,1:200);
TRIAL1 = FITNESS{num_gen+1};
% figure
% for i = 1:num_gen+1
%     plot(FITNESS{i}(:,2)*(10^3),FITNESS{i}(:,1),'.')
%     hold on
%     pause(0.01)
% end
% plot(TRIAL1(:,2)*(10^3),TRIAL1(:,1),'k')
load('maxFbkh_maxdeltalm_5000gen_rect_lateral_mar19_6x12x1_0.3E_0.5P_0.3M_0.1667h0_trial2_C01_260.5543e3.mat')
POPULATION_TRIAL2 = POPULATION(num_gen+1,1:200);
TRIAL2 = FITNESS{num_gen+1};
load('maxFbkh_maxdeltalm_5000gen_rect_lateral_mar19_6x12x1_0.3E_0.5P_0.3M_0.1667h0_trial3_C01_260.5543e3.mat')
POPULATION_TRIAL3 = POPULATION(num_gen+1,1:200);
TRIAL3 = FITNESS{num_gen+1};
load('maxFbkh_maxdeltalm_5000gen_rect_lateral_mar19_6x12x1_0.3E_0.5P_0.3M_0.1667h0_trial4_C01_260.5543e3.mat')
POPULATION_TRIAL4 = POPULATION(num_gen+1,1:200);
TRIAL4 = FITNESS{num_gen+1};
load('maxFbkh_maxdeltalm_5000gen_rect_lateral_mar19_6x12x1_0.3E_0.5P_0.3M_0.1667h0_trial5_C01_260.5543e3.mat')
POPULATION_TRIAL5 = POPULATION(num_gen+1,1:200);
TRIAL5 = FITNESS{num_gen+1};

POPULATION_LATERAL = cell(1,1000);
POPULATION_LATERAL(1:200) = POPULATION_TRIAL1;
POPULATION_LATERAL(201:400) = POPULATION_TRIAL2;
POPULATION_LATERAL(401:600) = POPULATION_TRIAL3;
POPULATION_LATERAL(601:800) = POPULATION_TRIAL4;
POPULATION_LATERAL(801:1000) = POPULATION_TRIAL5;
TRIALS_LAT_X = [TRIAL1(:,2)*(10^3); TRIAL2(:,2)*(10^3); TRIAL3(:,2)*(10^3); TRIAL4(:,2)*(10^3); TRIAL5(:,2)*(10^3)];
TRIALS_LAT_Y = [TRIAL1(:,1); TRIAL2(:,1); TRIAL3(:,1); TRIAL4(:,1); TRIAL5(:,1)];
TRIAL = [TRIALS_LAT_X, TRIALS_LAT_Y];
[membership_lat, member_value]=find_pareto_frontier(-TRIAL);
plot(-member_value(:,1),-member_value(:,2),'g')
hold on
c = 0;
for i = 1:1000
    if membership_lat(i) == 1
        c = c + 1;
        PARETOPOP_LAT{c} = POPULATION_LATERAL{i};
        PARETO_LATX(c) = TRIALS_LAT_X(i);
        PARETO_LATY(c) = TRIALS_LAT_Y(i);
        PARETOPOP_LAT_beta0(c) = PARETOPOP_LAT{c}(1,2); % [deg]
        PARETOPOP_LAT_l0(c) = PARETOPOP_LAT{c}(1,3)*2.54; % [cm]
        PARETOPOP_LAT_n(c) = 2*length(PARETOPOP_LAT{c}(:,1)); 
        PARETOPOP_LAT_r0(c) = PARETOPOP_LAT{c}(1,1)*25.4; % [mm]
    end
end
% 
% figure
load('maxFbkh_maxdeltalm_5000gen_rect_central_mar19_6x12x1_0.3E_0.5P_0.3M_0.1667h0_trial1_C01_260.5543e3.mat')
POPULATION_TRIAL1c = POPULATION(num_gen+1,1:200);
TRIAL1c = FITNESS{num_gen+1};
% for i = 1:num_gen+1
%     plot(FITNESS{i}(:,2)*(10^3),FITNESS{i}(:,1),'.')
%     hold on
% end
% plot(TRIAL1c(:,2)*(10^3),TRIAL1c(:,1),'k')
load('maxFbkh_maxdeltalm_5000gen_rect_central_mar19_6x12x1_0.3E_0.5P_0.3M_0.1667h0_trial2_C01_260.5543e3.mat')
POPULATION_TRIAL2c = POPULATION(num_gen+1,1:200);
TRIAL2c = FITNESS{num_gen+1};
% for i = 1:num_gen+1
%     plot(FITNESS{i}(:,2)*(10^3),FITNESS{i}(:,1),'.')
%     hold on
%     pause(0.01)
% end
% plot(TRIAL2c(:,2)*(10^3),TRIAL2c(:,1),'k')
load('maxFbkh_maxdeltalm_5000gen_rect_central_mar19_6x12x1_0.3E_0.5P_0.3M_0.1667h0_trial3_C01_260.5543e3.mat')
POPULATION_TRIAL3c = POPULATION(num_gen+1,1:200);
TRIAL3c = FITNESS{num_gen+1};
load('maxFbkh_maxdeltalm_5000gen_rect_central_mar19_6x12x1_0.3E_0.5P_0.3M_0.1667h0_trial4_C01_260.5543e3.mat')
POPULATION_TRIAL4c = POPULATION(num_gen+1,1:200);
TRIAL4c = FITNESS{num_gen+1};
load('maxFbkh_maxdeltalm_5000gen_rect_central_mar19_6x12x1_0.3E_0.5P_0.3M_0.1667h0_trial5_C01_260.5543e3.mat')
POPULATION_TRIAL5c = POPULATION(num_gen+1,1:200);
TRIAL5c = FITNESS{num_gen+1};

POPULATION_CENTRAL = cell(1,1000);
POPULATION_CENTRAL(1:200) = POPULATION_TRIAL1c;
POPULATION_CENTRAL(201:400) = POPULATION_TRIAL2c;
POPULATION_CENTRAL(401:600) = POPULATION_TRIAL3c;
POPULATION_CENTRAL(601:800) = POPULATION_TRIAL4c;
POPULATION_CENTRAL(801:1000) = POPULATION_TRIAL5c;
TRIALS_CEN_X = [TRIAL1c(:,2)*(10^3); TRIAL2c(:,2)*(10^3); TRIAL3c(:,2)*(10^3); TRIAL4c(:,2)*(10^3); TRIAL5c(:,2)*(10^3)];
TRIALS_CEN_Y = [TRIAL1c(:,1); TRIAL2c(:,1); TRIAL3c(:,1); TRIAL4c(:,1); TRIAL5c(:,1)];
TRIAL_C = [TRIALS_CEN_X, TRIALS_CEN_Y];
[membership_cen, member_value_central]=find_pareto_frontier(-TRIAL_C);
plot(-member_value_central(:,1),-member_value_central(:,2),'c')
c = 0;
for i = 1:1000
    if membership_cen(i) == 1
        c = c + 1;
        PARETOPOP_CEN{c} = POPULATION_CENTRAL{i};
        PARETO_CENX(c) = TRIALS_CEN_X(i);
        PARETO_CENY(c) = TRIALS_CEN_Y(i);
        PARETOPOP_CEN_beta0(c) = PARETOPOP_CEN{c}(1,2);
        PARETOPOP_CEN_l0(c) = PARETOPOP_CEN{c}(1,3)*2.54; % [cm]
        PARETOPOP_CEN_n(c) = 2*length(PARETOPOP_CEN{c}(:,1));
        PARETOPOP_CEN_r0(c) = PARETOPOP_CEN{c}(1,1)*25.4; % [mm]
    end
end
% close all
% figure
% plot(-member_value(:,1),-member_value(:,2),'g','LineWidth',1.5)
% hold on
% plot(-member_value_central(:,1),-member_value_central(:,2),'c','LineWidth',1.5)
% hold on
% plot(-member_value(end-8,1),-member_value(end-8,2),'ko')
% hold on
% plot(-member_value(end-9,1),-member_value(end-9,2),'k*')
% hold on
% plot(-member_value_central(end-29,1),-member_value_central(end-29,2),'ko')
% hold on
% plot(-member_value_central(end-30,1),-member_value_central(end-30,2),'k*')
% hold on
P2 = polyfit(-member_value(end-9:end-8,1),-member_value(end-9:end-8,2),1);
latx = -member_value(end-8,1):-member_value(end-9,1);
latfit = @(x) P2(1)*x+P2(2);
P = polyfit([-member_value_central(end-29,1) -member_value_central(end-30,1)],[-member_value_central(end-29,2) -member_value_central(end-30,2)],1);
cenx = -member_value_central(end-29,1):-member_value_central(end-30,1);
cenfit = @(x) P(1)*x+P(2);
xintersect = fzero(@(x) latfit(x)-cenfit(x),43);
% fplot(latfit,[-member_value(end-8,1) xintersect],'k-')
% hold on
% fplot(cenfit,[xintersect -member_value_central(end-30,1)],'k-')
% fsolve(@(x) latx(x)-cenfit(x),50)
% 
n = W/D; % [--]
P = 50*6894.76; % [Pa] 
F_piston = (W*0.0254)*(D*0.0254)*P;

% PARETO FRONTIER PLOT TO IDENTIFY DESIGNS ONLY
% figure
% plot(TRIAL1(:,2)*(10^3),TRIAL1(:,1),'r^')
% hold on
% plot(TRIALS_LAT_X,TRIALS_LAT_Y,'r^')
% hold on
% plot(-member_value(:,1),-member_value(:,2),'g^-','LineWidth',1.5)
% hold on
% plot(TRIALS_CEN_X,TRIALS_CEN_Y,'k^')
% hold on
% plot(-member_value_central(:,1),-member_value_central(:,2),'b^-','LineWidth',1.5)
% hold on
% ylim([0 Inf])
% grid on
% xlabel('\Deltal_{m,free} (mm)')
% ylabel('F_{b,mr} (N)')
% title('6 x 12 x 1')
% figure
% plot(-member_value(:,1)/(L*25.4),-member_value(:,2)/F_piston,'g^-','LineWidth',1.5)
% hold on
% plot(-member_value_central(:,1)/(L*25.4),-member_value_central(:,2)/F_piston,'b^-','LineWidth',1.5)
% hold on
% grid on
% xlabel('\Deltal_{m,free}/L')
% ylabel('F_{b,mr}/WDP ')
% title('6 x 12 x 1')


%% DOUBLE AXES PLOT OF PARETO FRONTIER WITH NON-DIMENSIONAL AXES
close all
figure
hAX = axes;
plot([-member_value(:,1); 0],[-member_value(:,2); max(-member_value(:,2))],'g-','LineWidth',3)
hold on
plot([-member_value_central(:,1); 0],[-member_value_central(:,2); max(-member_value_central(:,2))],'c-','LineWidth',3)
hold on
plot([-member_value(end,1); 0],[-member_value(end,2); -member_value(end,2)],'k-','LineWidth',1.5)
hold on
plot([-member_value(end-8:end,1); 0],[-member_value(end-8:end,2); max(-member_value(end-8:end,2))],'k-','LineWidth',1.5)
hold on
plot([-member_value(end-8,1) xintersect],[-member_value(end-8,2) cenfit(xintersect)],'k-','LineWidth',1.5)
hold on
plot([xintersect -member_value_central(end-30,1)],[cenfit(xintersect) -member_value_central(end-30,2)],'k-','LineWidth',1.5)
hold on
plot(-member_value_central(1:end-30,1),-member_value_central(1:end-30,2),'k-','LineWidth',1.5)
hold on
plot([-member_value_central(1,1) -member_value_central(1,1)],[-member_value_central(1,2) 0],'k-','LineWidth',1.5)
grid on
xlim([0 Inf])
ylim([0 Inf])
title('6 x 12 x 1')
pos = get(hAX,'position');
pos_left = pos(1);
pos(1) = pos_left*2; 
pos(3) = pos(3)-pos(1);
set(hAX(1),'position',pos)
pos(1) = pos(1)-pos_left;
pos(3) = 0.01;
pos(2) = pos(2)+pos(2);
pos(4) = pos(4)-0.5*pos(2);
hAX(2) = axes('XTickLabel',{},'XTick',[],'YTickMode','manual','YLimMode','manual','YLim',[0 max(max(-member_value(:,2))/F_piston,max(-member_value_central(:,2))/F_piston)],'TickLength',[0.01 pos(3)],'YTick',[0:0.1:1.5],'XTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1','1.1','1.2','1.3','1.4','1.5'},'position',pos,'color','none');
% hAX(2) = axes('XTickLabel',{},'XTick',[],'YTickMode','manual','YLimMode','manual','YLim',[0 max(max(-member_value(:,2))/F_piston,max(-member_value_central(:,2))/F_piston)],'TickLength',[0.01 pos(3)],'YTick',[0:0.5:1.5],'XTickLabel',{'0','0.5','1','1.5'},'position',pos,'color','none');
set(hAX(2),'position',pos)
pos = get(hAX(1),'position');
pos1 = pos(2);
pos(2) = pos(2)+pos1; pos(4) = pos(4)-pos1;
set(hAX(1),'position',pos)
pos(2) = pos1-0.01; pos(4) = 0.001;
hAX(3) = axes('XTickMode','manual','XLimMode','manual','XLim',[0 max(max(-member_value(:,1))/(L*25.4),max(-member_value_central(:,1))/(L*25.4))],'TickLength',[0.01 pos(4)],'XTick',[0:0.1:0.8],'XTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'},'position',pos,'color','none');
grid on
xlabel(hAX(1),'\Deltal_{m,free} (mm)','FontSize',9)
ylabel(hAX(1),'F_{b,MR} (N)','FontSize',9)
ylabel(hAX(2),'^{F_{b,MR}}/{WDP}')
xlabel(hAX(3),'\Deltal_{m,free}/L')
%%
close all
% DIAGNOSTIC PLOTS
% PARAMETERS VS. BUNDLE CONTRACTION
figure
[sorted_PARETO_CENX,Index_PARETO_CENX] = sort(PARETO_CENX);
plot(sorted_PARETO_CENX(38:end),PARETO_CENY(Index_PARETO_CENX((38:end))),'o')
hold on
[sorted_PARETO_LATX,Index_PARETO_LATX] = sort(PARETO_LATX);
plot(sorted_PARETO_LATX(1:19),PARETO_LATY(Index_PARETO_LATX(1:19)),'o')
hold on

%%
close all
figure
hAX = axes;
plot(sorted_PARETO_LATX(1:19),PARETOPOP_LAT_beta0(Index_PARETO_LATX(1:19)),'bo')
hold on
plot(sorted_PARETO_LATX(1:19),PARETOPOP_LAT_l0(Index_PARETO_LATX(1:19))*2,'rs')
hold on
plot(sorted_PARETO_LATX(1:19),PARETOPOP_LAT_n(Index_PARETO_LATX(1:19))*5,'^','MarkerEdgeColor',[0 0.5 0])
hold on
xline(sorted_PARETO_LATX(19),'k--','LineWidth',1.5)
hold on
plot(sorted_PARETO_CENX(38:end),PARETOPOP_CEN_beta0(Index_PARETO_CENX(38:end)),'bo')
hold on
plot(sorted_PARETO_CENX(38:end),PARETOPOP_CEN_l0(Index_PARETO_CENX(38:end))*2,'rs')
hold on
plot(sorted_PARETO_CENX(38:end),PARETOPOP_CEN_n(Index_PARETO_CENX(38:end))*5,'^','MarkerEdgeColor',[0 0.5 0])
hold on
grid on
xlim([0 Inf])
ylim([0 Inf])
% title('6 x 12 x 1')
pos = get(hAX,'position');
pos_left = pos(1);
pos(1) = pos_left*3; 
pos(3) = pos(3)-pos(1);
set(hAX(1),'position',pos)
pos(1) = pos(1)-pos_left;
pos(3) = 0.01;
pos(2) = pos(2);
hAX(2) = axes('XTickLabel',{},'XTick',[],'YTickMode','manual','YLimMode','manual','TickLength',[0.01 pos(3)],'YTick',[0:(1/6):0.95],'YTickLabel',{'0','5','10','15','20','25'},'position',pos,'color','none','YColor','r');
set(hAX(2),'position',pos)
pos(1) = pos(1)-pos_left;
pos(3) = 0.01;
pos(2) = pos(2);
hAX(3) = axes('XTickLabel',{},'XTick',[],'YTickMode','manual','YLimMode','manual','TickLength',[0.01 pos(3)],'YTick',[0:(1/6):0.95],'YTickLabel',{'0','2','4','6','8','10'},'position',pos,'color','none','YColor',[0 0.5 0]);
xlabel(hAX(1),'\Deltal_{m,free} (mm)','FontSize',10,'FontWeight','bold')
ylabel(hAX(1),'\beta_0 (deg.)','FontSize',10,'Color','b','FontWeight','bold')
ylabel(hAX(2),'l_0 (cm)','FontSize',10,'FontWeight','bold')
ylabel(hAX(3),'n','FontSize',10,'FontWeight','bold')
%%
close all
% PARAMETERS VS. BUNDLE FORCE
figure
[sorted_PARETO_CENY,Index_PARETO_CENY] = sort(PARETO_CENY);
plot(PARETO_CENX(Index_PARETO_CENY(1:end-39)),sorted_PARETO_CENY(1:end-39),'o')
hold on
[sorted_PARETO_LATY,Index_PARETO_LATY] = sort(PARETO_LATY);
plot(PARETO_LATX(Index_PARETO_LATY(end-16:end)),sorted_PARETO_LATY(end-16:end),'o')
%%
close all
figure
hAX = axes;
plot(sorted_PARETO_CENY(1:end-39),PARETOPOP_CEN_beta0(Index_PARETO_CENY(1:end-39)),'bo')
hold on
plot(sorted_PARETO_CENY(1:end-39),PARETOPOP_CEN_l0(Index_PARETO_CENY(1:end-39))*2,'rs')
hold on
plot(sorted_PARETO_CENY(1:end-39),PARETOPOP_CEN_n(Index_PARETO_CENY(1:end-39))*5,'^','MarkerEdgeColor',[0 0.5 0])
hold on
xline(sorted_PARETO_CENY(end-39),'k--','LineWidth',1.5)
hold on
plot(sorted_PARETO_LATY(end-16:end),PARETOPOP_LAT_beta0(Index_PARETO_LATY(end-16:end)),'bo')
hold on
plot(sorted_PARETO_LATY(end-16:end),PARETOPOP_LAT_l0(Index_PARETO_LATY(end-16:end))*2,'rs')
hold on
plot(sorted_PARETO_LATY(end-16:end),PARETOPOP_LAT_n(Index_PARETO_LATY(end-16:end))*5,'^','MarkerEdgeColor',[0 0.5 0])
hold on
grid on
xlim([0 Inf])
ylim([0 Inf])
% title('6 x 12 x 1')
pos = get(hAX,'position');
pos_left = pos(1);
pos(1) = pos_left*3; 
pos(3) = pos(3)-pos(1);
set(hAX(1),'position',pos)
pos(1) = pos(1)-pos_left;
pos(3) = 0.01;
pos(2) = pos(2);
hAX(2) = axes('XTickLabel',{},'XTick',[],'YTickMode','manual','YLimMode','manual','TickLength',[0.01 pos(3)],'YTick',[0:(1/6):0.95],'YTickLabel',{'0','5','10','15','20','25'},'position',pos,'color','none','YColor','r');
set(hAX(2),'position',pos)
pos(1) = pos(1)-pos_left;
pos(3) = 0.01;
pos(2) = pos(2);
hAX(3) = axes('XTickLabel',{},'XTick',[],'YTickMode','manual','YLimMode','manual','TickLength',[0.01 pos(3)],'YTick',[0:(1/6):0.95],'YTickLabel',{'0','2','4','6','8','10'},'position',pos,'color','none','YColor',[0 0.5 0]);
xlabel(hAX(1),'F_{b,mr} (N)','FontSize',10,'FontWeight','bold')
ylabel(hAX(1),'\beta_0 (deg.)','FontSize',10,'Color','b','FontWeight','bold')
ylabel(hAX(2),'l_0 (cm)','FontSize',10,'FontWeight','bold')
ylabel(hAX(3),'n','FontSize',10,'FontWeight','bold')