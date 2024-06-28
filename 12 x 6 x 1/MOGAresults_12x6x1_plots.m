%% 12 x 6 x 1 
% LATERAL + CENTRAL PARETO FRONTIER 
clc
clear all
close all

load('maxFbkh_maxdeltalm_5000gen_rect_lateral_mar19_12x6x1_0.3E_0.5P_0.3M_0.1667h0_trial1_C01_260.5543e3.mat')
POPULATION_TRIAL1 = POPULATION(num_gen+1,1:200);
TRIAL1 = FITNESS{1:num_gen+1};
load('maxFbkh_maxdeltalm_5000gen_rect_lateral_mar19_12x6x1_0.3E_0.5P_0.3M_0.1667h0_trial2_C01_260.5543e3.mat')
POPULATION_TRIAL2 = POPULATION(num_gen+1,1:200);
TRIAL2 = FITNESS{1:num_gen+1};
load('maxFbkh_maxdeltalm_5000gen_rect_lateral_mar19_12x6x1_0.3E_0.5P_0.3M_0.1667h0_trial3_C01_260.5543e3.mat')
POPULATION_TRIAL3 = POPULATION(num_gen+1,1:200);
TRIAL3 = FITNESS{1:num_gen+1};
load('maxFbkh_maxdeltalm_5000gen_rect_lateral_mar19_12x6x1_0.3E_0.5P_0.3M_0.1667h0_trial4_C01_260.5543e3.mat')
POPULATION_TRIAL4 = POPULATION(num_gen+1,1:200);
TRIAL4 = FITNESS{1:num_gen+1};
load('maxFbkh_maxdeltalm_5000gen_rect_lateral_mar19_12x6x1_0.3E_0.5P_0.3M_0.1667h0_trial5_C01_260.5543e3.mat')
POPULATION_TRIAL5 = POPULATION(num_gen+1,1:200);
TRIAL5 = FITNESS{1:num_gen+1};

load('maxFbkh_maxdeltalm_5000gen_rect_central_mar19_12x6x1_0.3E_0.5P_0.3M_0.1667h0_trial1_C01_260.5543e3.mat')
POPULATION_TRIAL1c = POPULATION(num_gen+1,1:200);
TRIAL1c = FITNESS{1:num_gen+1};
load('maxFbkh_maxdeltalm_5000gen_rect_central_mar19_12x6x1_0.3E_0.5P_0.3M_0.1667h0_trial2_C01_260.5543e3.mat')
POPULATION_TRIAL2c = POPULATION(num_gen+1,1:200);
TRIAL2c = FITNESS{1:num_gen+1};
load('maxFbkh_maxdeltalm_5000gen_rect_central_mar19_12x6x1_0.3E_0.5P_0.3M_0.1667h0_trial3_C01_260.5543e3.mat')
POPULATION_TRIAL3c = POPULATION(num_gen+1,1:200);
TRIAL3c = FITNESS{1:num_gen+1};
load('maxFbkh_maxdeltalm_5000gen_rect_central_mar19_12x6x1_0.3E_0.5P_0.3M_0.1667h0_trial4_C01_260.5543e3.mat')
POPULATION_TRIAL4c = POPULATION(num_gen+1,1:200);
TRIAL4c = FITNESS{1:num_gen+1};
load('maxFbkh_maxdeltalm_5000gen_rect_central_mar19_12x6x1_0.3E_0.5P_0.3M_0.1667h0_trial5_C01_260.5543e3.mat')
POPULATION_TRIAL5c = POPULATION(num_gen+1,1:200);
TRIAL5c = FITNESS{1:num_gen+1};

% TRIALS_LAT_X = [TRIAL1(:,2); TRIAL2(:,2); TRIAL3(:,2); TRIAL4(:,2); TRIAL5(:,2)]*((10^2)^3);
TRIALS_LAT_X = [TRIAL1(:,2)*(10^3); TRIAL2(:,2)*(10^3); TRIAL3(:,2)*(10^3); TRIAL4(:,2)*(10^3); TRIAL5(:,2)*(10^3)];
TRIALS_LAT_Y = [TRIAL1(:,1); TRIAL2(:,1); TRIAL3(:,1); TRIAL4(:,1); TRIAL5(:,1)];
TRIAL = [TRIALS_LAT_X, TRIALS_LAT_Y];
[membership, member_value]=find_pareto_frontier(-TRIAL);

% TRIALS_CEN_X = [TRIAL1c(:,2); TRIAL2c(:,2); TRIAL3c(:,2); TRIAL4c(:,2); TRIAL5c(:,2)]*((10^2)^3);
TRIALS_CEN_X = [TRIAL1c(:,2)*(10^3); TRIAL2c(:,2)*(10^3); TRIAL3c(:,2)*(10^3); TRIAL4c(:,2)*(10^3); TRIAL5c(:,2)*(10^3)];
TRIALS_CEN_Y = [TRIAL1c(:,1); TRIAL2c(:,1); TRIAL3c(:,1); TRIAL4c(:,1); TRIAL5c(:,1)];
TRIAL_C = [TRIALS_CEN_X, TRIALS_CEN_Y];
[membership, member_value_central]=find_pareto_frontier(-TRIAL_C);

plot([max(-member_value(:,1)); -member_value(:,1); 0;],[0; -member_value(:,2); max(-member_value(:,2));],'g-','LineWidth',2)
hold on
plot([max(-member_value_central(:,1)); -member_value_central(:,1); 0;],[0; -member_value_central(:,2); max(-member_value_central(:,2));],'c-','LineWidth',2)
hold on

xlim([0 Inf])
grid on
xlabel('\Deltal_{m,free} (mm)')
ylabel('F_{b,MR} (N)')
title('Rectangular LxWxD = 12"x6"x1"')
%% LATERAL + CENTRAL PARETO FRONTIER (w/MARKERS FOR DESIGNS ALONG PARETO FRONTIER)
clc
clear all
close all

load('maxFbkh_maxdeltalm_5000gen_rect_lateral_mar19_12x6x1_0.3E_0.5P_0.3M_0.1667h0_trial1_C01_260.5543e3.mat')
POPULATION_TRIAL1 = POPULATION(num_gen+1,1:200);
TRIAL1 = FITNESS{num_gen+1};
load('maxFbkh_maxdeltalm_5000gen_rect_lateral_mar19_12x6x1_0.3E_0.5P_0.3M_0.1667h0_trial2_C01_260.5543e3.mat')
POPULATION_TRIAL2 = POPULATION(num_gen+1,1:200);
TRIAL2 = FITNESS{num_gen+1};
load('maxFbkh_maxdeltalm_5000gen_rect_lateral_mar19_12x6x1_0.3E_0.5P_0.3M_0.1667h0_trial3_C01_260.5543e3.mat')
POPULATION_TRIAL3 = POPULATION(num_gen+1,1:200);
TRIAL3 = FITNESS{num_gen+1};
load('maxFbkh_maxdeltalm_5000gen_rect_lateral_mar19_12x6x1_0.3E_0.5P_0.3M_0.1667h0_trial4_C01_260.5543e3.mat')
POPULATION_TRIAL4 = POPULATION(num_gen+1,1:200);
TRIAL4 = FITNESS{num_gen+1};
load('maxFbkh_maxdeltalm_5000gen_rect_lateral_mar19_12x6x1_0.3E_0.5P_0.3M_0.1667h0_trial5_C01_260.5543e3.mat')
POPULATION_TRIAL5 = POPULATION(num_gen+1,1:200);
TRIAL5 = FITNESS{num_gen+1};

load('maxFbkh_maxdeltalm_5000gen_rect_central_mar19_12x6x1_0.3E_0.5P_0.3M_0.1667h0_trial1_C01_260.5543e3.mat')
POPULATION_TRIAL1c = POPULATION(num_gen+1,1:200);
TRIAL1c = FITNESS{num_gen+1};
load('maxFbkh_maxdeltalm_5000gen_rect_central_mar19_12x6x1_0.3E_0.5P_0.3M_0.1667h0_trial2_C01_260.5543e3.mat')
POPULATION_TRIAL2c = POPULATION(num_gen+1,1:200);
TRIAL2c = FITNESS{num_gen+1};
load('maxFbkh_maxdeltalm_5000gen_rect_central_mar19_12x6x1_0.3E_0.5P_0.3M_0.1667h0_trial3_C01_260.5543e3.mat')
POPULATION_TRIAL3c = POPULATION(num_gen+1,1:200);
TRIAL3c = FITNESS{num_gen+1};
load('maxFbkh_maxdeltalm_5000gen_rect_central_mar19_12x6x1_0.3E_0.5P_0.3M_0.1667h0_trial4_C01_260.5543e3.mat')
POPULATION_TRIAL4c = POPULATION(num_gen+1,1:200);
TRIAL4c = FITNESS{num_gen+1};
load('maxFbkh_maxdeltalm_5000gen_rect_central_mar19_12x6x1_0.3E_0.5P_0.3M_0.1667h0_trial5_C01_260.5543e3.mat')
POPULATION_TRIAL5c = POPULATION(num_gen+1,1:200);
TRIAL5c = FITNESS{num_gen+1};

figure
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
plot(-member_value(:,1),-member_value(:,2),'g-','LineWidth',1.5)
hold on
c = 0;
for i = 1:1000
    if membership_lat(i) == 1
        c = c + 1;
        PARETOPOP_LAT{c} = POPULATION_LATERAL{i};
        PARETO_LATX(c) = TRIALS_LAT_X(i);
        PARETO_LATY(c) = TRIALS_LAT_Y(i);
        PARETOPOP_LAT_beta0(c) = PARETOPOP_LAT{c}(1,2);
        PARETOPOP_LAT_l0(c) = PARETOPOP_LAT{c}(1,3)*2.54;
        PARETOPOP_LAT_n(c) = 2*length(PARETOPOP_LAT{c}(:,1));
        PARETOPOP_LAT_r0(c) = PARETOPOP_LAT{c}(1,1)*25.4; 
        plot(PARETO_LATX(c),PARETO_LATY(c),'g^','LineWidth',1.5)
        hold on
        % plot(PARETO_LATX(c),PARETOPOP_LAT{c}(1,3)*2.54,'ks','MarkerFaceColor',[0.4660 0.6740 0.1880])
        % plot(PARETO_LATX(c),PARETOPOP_LAT{c}(1,2),'ks','MarkerFaceColor',[0.4660 0.6740 0.1880])
        % plot(PARETO_LATY(c),PARETOPOP_LAT{c}(1,1)*25.4,'ks','MarkerFaceColor',[0.4660 0.6740 0.1880])
        % plot(PARETO_LATY(c),PARETOPOP_LAT{c}(1,2),'ks','MarkerFaceColor',[0.4660 0.6740 0.1880])
        % plot(PARETO_LATY(c),2*length(PARETOPOP_LAT{c}(:,1)),'ks','MarkerFaceColor',[0.4660 0.6740 0.1880])
        % hold on
    end
end
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
plot(-member_value_central(:,1),-member_value_central(:,2),'b-','LineWidth',1.5)
hold on

c = 0;
for i = 1:1000
    if membership_cen(i) == 1
        c = c + 1;
        PARETOPOP_CEN{c} = POPULATION_CENTRAL{i};
        PARETO_CENX(c) = TRIALS_CEN_X(i);
        PARETO_CENY(c) = TRIALS_CEN_Y(i);
        PARETOPOP_CEN_beta0(c) = PARETOPOP_CEN{c}(1,2);
        PARETOPOP_CEN_l0(c) = PARETOPOP_CEN{c}(1,3)*2.54;
        PARETOPOP_CEN_n(c) = 2*length(PARETOPOP_CEN{c}(:,1));
        PARETOPOP_CEN_r0(c) = PARETOPOP_CEN{c}(1,1)*25.4; 
        plot(PARETO_CENX(c),PARETO_CENY(c),'b^','LineWidth',1.5)
        hold on
        % plot(PARETO_CENX(c),PARETOPOP_CEN{c}(1,3)*2.54,'ks','MarkerFaceColor','b')
        % plot(PARETO_CENX(c),PARETOPOP_CEN{c}(1,2),'ks','MarkerFaceColor','b')
        % plot(PARETO_CENY(c),PARETOPOP_CEN{c}(1,1)*25.4,'ks','MarkerFaceColor','b')
        % plot(PARETO_CENY(c),PARETOPOP_CEN{c}(1,2),'ks','MarkerFaceColor','b')
        % plot(PARETO_CENY(c),2*length(PARETOPOP_CEN{c}(:,1)),'ks','MarkerFaceColor','b')
        % hold on
    end
end
grid on
xlim([0 Inf])
grid on
xlabel('\Deltal_{m,free} (mm)')
ylabel('F_{b,MR} (N)')
title('Rectangular LxWxD = 12"x6"x1"')

%% DOUBLE AXES (F vs. deltalmfree, F/WDP vs deltalmfree/L)
clc
close all
figure
hAX = axes;
plot([0 -member_value(end,1)],[-member_value(end,2) -member_value(end,2)],'g-','LineWidth',3)
hold on
plot(-member_value(:,1),-member_value(:,2),'g-','LineWidth',3)
hold on
plot([0 -member_value_central(end,1)],[-member_value_central(end,2) -member_value_central(end,2)],'c-','Linewidth',3)
hold on
plot(-member_value_central(:,1),-member_value_central(:,2),'c-','LineWidth',3)
hold on
plot([0 -member_value_central(end,1)],[-member_value_central(end,2) -member_value_central(end,2)],'k-','Linewidth',1.5)
hold on
plot(-member_value_central(17:end,1),-member_value_central(17:end,2),'k-','LineWidth',1.5)
hold on
plot(-member_value(:,1),-member_value(:,2),'k-','LineWidth',1.5)
hold on
P1 = polyfit(-member_value_central(16:17,1),-member_value_central(16:17,2),1);
cenx = -member_value_central(17,1):-member_value_central(16,1);
cenfit = @(x) P1(1)*x+P1(2);
P2 = polyfit([0 -member_value(end,1)],[-member_value(end,2) -member_value(end,2)],1);
latx = 0:-member_value(end,1);
latfit = @(x) P2(1)*x+P2(2);
xintersect = fzero(@(x) latfit(x)-cenfit(x),61);
fplot(cenfit,[-member_value_central(17,1) xintersect],'k-','LineWidth',1.5)
hold on
fplot(latfit,[xintersect -member_value(end,1)],'k-','LineWidth',1.5)
hold on
grid on
xlim([0 Inf])
ylim([0 Inf])
title('12 x 6 x 1')
pos = get(hAX,'position');
pos_left = pos(1);
pos(1) = pos_left*2; 
pos(3) = pos(3)-pos(1);
set(hAX(1),'position',pos)
% set(hAX(1),'xticklabel',[], 'yticklabel', [])
P = 50*6894.76; % [Pa] 
F_piston = (W*0.0254)*(D*0.0254)*P;
pos(1) = pos(1)-pos_left;
pos(3) = 0.01;
pos(2) = pos(2)+pos(2);
pos(4) = pos(4)-0.5*pos(2);
hAX(2) = axes('XTickLabel',{},'XTick',[],'YTickMode','manual','YLimMode','manual','YLim',[0 max(max(-member_value(:,2))/F_piston,max(-member_value_central(:,2))/F_piston)],'TickLength',[0.01 pos(3)],'YTick',[0:0.5:4.5],'XTickLabel',{'0','0.5','1','1.5','2','2.5','3','3.5','4','4.5'},'position',pos,'color','none');
set(hAX(2),'position',pos)
pos = get(hAX(1),'position');
pos1 = pos(2);
pos(2) = pos(2)+pos1; pos(4) = pos(4)-pos1;
set(hAX(1),'position',pos)
pos(2) = pos1-0.01; pos(4) = 0.001;
hAX(3) = axes('XTickMode','manual','XLimMode','manual','XLim',[0 max(max(-member_value(:,1))/(L*25.4),max(-member_value_central(:,1))/(L*25.4))],'TickLength',[0.01 pos(4)],'XTick',[0:0.05:0.35],'XTickLabel',{'0','0.05','0.1','0.15','0.2','0.25','0.3','0.35'},'position',pos,'color','none');
grid on
xlabel(hAX(1),'\Deltal_{m,free} (mm)','FontSize',9)
ylabel(hAX(1),'F_{b,MR} (N)','FontSize',9)
ylabel(hAX(2),'^{F_{b,MR}}/{WDP}')
xlabel(hAX(3),'\Deltal_{m,free}/L')
% xlabel('\Deltal_{free} (mm)')
% ylabel('F_{b,MR} (N)')
% title('12 x 6 x 1')
%% DIAGNOSTIC PLOTS
close all
% PARAMETERS VS. BUNDLE CONTRACTION
% figure
[sorted_PARETO_CENX,Index_PARETO_CENX] = sort(PARETO_CENX);
plot(sorted_PARETO_CENX(1:end-32),PARETO_CENY(Index_PARETO_CENX((1:end-32))),'o')
hold on
[sorted_PARETO_LATX,Index_PARETO_LATX] = sort(PARETO_LATX);
plot(sorted_PARETO_LATX(1:end),PARETO_LATY(Index_PARETO_LATX(1:end)),'o')
hold on
%%
close all
figure
hAX = axes;
plot(sorted_PARETO_CENX(1:end-32),PARETOPOP_CEN_beta0(Index_PARETO_CENX(1:end-32)),'bo')
hold on
plot(sorted_PARETO_CENX(1:end-32),PARETOPOP_CEN_l0(Index_PARETO_CENX(1:end-32)),'rs')
hold on
plot(sorted_PARETO_CENX(1:end-32),PARETOPOP_CEN_n(Index_PARETO_CENX(1:end-32))*2,'^','MarkerEdgeColor',[0, 0.5, 0])
hold on
xline(sorted_PARETO_CENX(end-32),'k--','LineWidth',1.5)
hold on
plot(sorted_PARETO_LATX(1:end),PARETOPOP_LAT_beta0(Index_PARETO_LATX(1:end)),'bo')
hold on
plot(sorted_PARETO_LATX(1:end),PARETOPOP_LAT_l0(Index_PARETO_LATX(1:end)),'rs')
hold on
plot(sorted_PARETO_LATX(1:end),PARETOPOP_LAT_n(Index_PARETO_LATX(1:end))*2,'^','MarkerEdgeColor',[0, 0.5, 0])
hold on
grid on
xlim([0 Inf])
ylim([0 60])
title('12 x 6 x 1')
pos = get(hAX,'position');
pos_left = pos(1);
pos(1) = pos_left*3; 
pos(3) = pos(3)-pos(1);
set(hAX(1),'position',pos)
pos(1) = pos(1)-pos_left;
pos(3) = 0.01;
pos(2) = pos(2);
hAX(2) = axes('XTickLabel',{},'XTick',[],'YTickMode','manual','YLimMode','manual','TickLength',[0.01 pos(3)],'YTick',[0:(1/6):0.95],'YTickLabel',{'0','10','20','30','40','50'},'position',pos,'color','none','YColor','r');
set(hAX(2),'position',pos)
pos(1) = pos(1)-pos_left;
pos(3) = 0.01;
pos(2) = pos(2);
hAX(3) = axes('XTickLabel',{},'XTick',[],'YTickMode','manual','YLimMode','manual','TickLength',[0.01 pos(3)],'YTick',[0:(1/15):0.95],'YTickLabel',{'0','2','4','6','8','10','12','14','16','18','20','22','24','26','28'},'position',pos,'color','none','YColor',[0, 0.5, 0]);

% hAX(3) = axes('XTickLabel',{},'XTick',[],'YTickMode','manual','YLimMode','manual','TickLength',[0.01 pos(3)],'YTick',[0:(1/6):0.95],'YTickLabel',{'0','5','10','15','20','25'},'position',pos,'color','none','YColor',[0, 0.5, 0]);
xlabel(hAX(1),'\Deltal_{m,free} (mm)','FontSize',10,'FontWeight','bold')
ylabel(hAX(1),'\beta_0 (deg.)','FontSize',10,'FontWeight','bold','Color','b')
ylabel(hAX(2),'l_0 (cm)','FontSize',10,'FontWeight','bold')
ylabel(hAX(3),'n','FontSize',10,'FontWeight','bold')
%%
close all
% PARAMETERS VS. BUNDLE FORCE
figure
[sorted_PARETO_CENY,Index_PARETO_CENY] = sort(PARETO_CENY);
plot(PARETO_CENX(Index_PARETO_CENY(35:end)),sorted_PARETO_CENY(35:end),'o')
hold on
[sorted_PARETO_LATY,Index_PARETO_LATY] = sort(PARETO_LATY);
plot(PARETO_LATX(Index_PARETO_LATY(1:end)),sorted_PARETO_LATY(1:end),'o')

%%
close all
figure
hAX = axes;
plot(sorted_PARETO_CENY(35:end),PARETOPOP_CEN_beta0(Index_PARETO_CENY(35:end)),'bo')
hold on
plot(sorted_PARETO_CENY(35:end),PARETOPOP_CEN_l0(Index_PARETO_CENY(35:end)),'rs')
hold on
plot(sorted_PARETO_CENY(35:end),PARETOPOP_CEN_n(Index_PARETO_CENY(35:end))*2,'^','MarkerEdgeColor',[0, 0.5, 0])
hold on
plot(sorted_PARETO_LATY(1:end),PARETOPOP_LAT_beta0(Index_PARETO_LATY(1:end)),'bo')
hold on
plot(sorted_PARETO_LATY(1:end),PARETOPOP_LAT_l0(Index_PARETO_LATY(1:end)),'rs')
hold on
plot(sorted_PARETO_LATY(1:end),PARETOPOP_LAT_n(Index_PARETO_LATY(1:end))*2,'^','MarkerEdgeColor',[0, 0.5, 0])
hold on
xline(sorted_PARETO_LATY(end),'k--','LineWidth',1.5)
hold on
% plot([1005.5 2001.1 3016.6 4022.1],ones(1,4)*PARETOPOP_CEN_beta0(Index_PARETO_CENY(1)),'ko')
% hold on
% plot([1005.5 2001.1 3016.6 4022.1],ones(1,4)*PARETOPOP_CEN_l0(Index_PARETO_CENY(1))*4,'rs')
% hold on
% plot([1005.5 2001.1 3016.6 4022.1],ones(1,4)*PARETOPOP_CEN_n(Index_PARETO_CENY(1)),'^','MarkerEdgeColor',[0.4660 0.6740 0.1880])
% hold on
grid on
xlim([0 Inf])
ylim([0 Inf])
title('12 x 6 x 1')
pos = get(hAX,'position');
pos_left = pos(1);
pos(1) = pos_left*3; 
pos(3) = pos(3)-pos(1);
set(hAX(1),'position',pos)
pos(1) = pos(1)-pos_left;
pos(3) = 0.01;
pos(2) = pos(2);
hAX(2) = axes('XTickLabel',{},'XTick',[],'YTickMode','manual','YLimMode','manual','TickLength',[0.01 pos(3)],'YTick',[0:(1/6):0.95],'YTickLabel',{'0','10','20','30','40','50'},'position',pos,'color','none','YColor','r');
set(hAX(2),'position',pos)
pos(1) = pos(1)-pos_left;
pos(3) = 0.01;
pos(2) = pos(2);
hAX(3) = axes('XTickLabel',{},'XTick',[],'YTickMode','manual','YLimMode','manual','TickLength',[0.01 pos(3)],'YTick',[0:(1/15):0.95],'YTickLabel',{'0','2','4','6','8','10','12','14','16','18','20','22','24','26','28'},'position',pos,'color','none','YColor',[0, 0.5, 0]);
% hAX(3) = axes('XTickLabel',{},'XTick',[],'YTickMode','manual','YLimMode','manual','TickLength',[0.01 pos(3)],'YTick',[0:(1/6):0.95],'YTickLabel',{'0','5','10','15','20','25'},'position',pos,'color','none','YColor',[0, 0.5, 0]);
xlabel(hAX(1),'F_{b,MR} (N)','FontSize',10,'FontWeight','bold')
ylabel(hAX(1),'\beta_0 (deg.)','FontSize',10,'FontWeight','bold','Color','b')
ylabel(hAX(2),'l_0 (cm)','FontSize',10,'FontWeight','bold')
ylabel(hAX(3),'n','FontSize',10,'FontWeight','bold')
