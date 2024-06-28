clc
clear all
close all 
%% DESIGN 3
load('June_3.5in_50psi_easier2plot.mat')
plot(DATA(:,2)*25.4,DATA(:,1)*1000*4.448) % Displacement (mm), Force (lbf)
hold on

P = 50*6894.76; % [Pa]
h0 = (1/16)/(3/8);
n = 2;
l0 = 3.5*0.0254;
beta0 = 51;
r0 = 0.3749*0.0254;
alpha0 = 32;

beta = @(alpha) real(asind((sind(beta0)*cosd(alpha0))/cosd(alpha)));
l = @(alpha) l0*(cosd(alpha)/cosd(alpha0));
C10 = 55.4298e3;
C01 = 260.5543e3;
t0 = h0*(r0);
B = l0/cosd(alpha0);     % [m] - length of fiber cord around FAM
N = (B*sind(alpha0))/(2*pi*r0); % [--] - number of fiber turns around FAM
F_gaylord = @(lambda) (P*((3*((lambda*l0)^2))-(B^2)))/(4*(N^2)*pi); 
Vb = (pi*(r0^2)*l0) - (pi*((r0-t0)^2)*l0); % [m^3] - bladder volume
LAMBDA = @(alpha) cosd(alpha)/cosd(alpha0);
term1 = @(lambda) (4*(C10+C01))*((l0^2)*(-1+(lambda^4)));
term2 = @(lambda) (4.*(l0.^6)*(-1+lambda)*(lambda^2)*(1+lambda)*(C10+(C01*(lambda^2))))/(((-4*(N^2)*(pi^2)*(r0^2))+((l0^2)*(-1+(lambda^2))))^2);
term3 = @(lambda) ((4*(l0^4)*(C10+(C01*(lambda^4))))/((-4*(N^2)*(pi^2)*(r0^2))+((l0^2)*(-1+(lambda^2)))));
term4 = @(lambda) ((l0^4)*(lambda^4)*(C10+(C01*(-1+(2*(lambda.^2))))))/((N^2)*(pi^2)*(r0^2));
F = @(lambda) F_gaylord(lambda) + Vb*((1./(2*(l0^3)*(lambda^3))))*(term1(lambda)+term2(lambda)-term3(lambda)-term4(lambda));
F = @(alpha) n*F(LAMBDA(alpha)).*cosd(beta(alpha)); % For each FAM pair in the bundle
alpha_range = alpha0:0.01:atand(sqrt(2));
for i = 1:length(alpha_range)
    model_deltalm2(i) = ((l0*cosd(beta0))-(l(alpha_range(i))*cosd(beta(alpha_range(i)))))/0.0254;
    model_F2(i) = F(alpha_range(i));
end
% alpha_c = fsolve(@(alpha) F(alpha),alpha0);
Fb = F(alpha0)*0.2248089431; % Blocked force (lbf)
alpha_free = fsolve(F,alpha0); % Braid angle at bundle free contraction (deg.)
% FAM rotation behavior changes due to inactive length
l_in = (10.5*25.4)/1000;
beta = @(alpha) asind((sind(beta0)*(l0+l_in))/(l(alpha)+l_in));
deltalm_free = (((l0+l_in)*cosd(beta0)) - ((l_in+l(alpha_free))*cosd(beta(alpha_free))))*1000; % Bundle free contraction [mm]
fplot(@(alpha) (((l0+l_in)*cosd(beta0)) - ((l_in+l(alpha))*cosd(beta(alpha))))*1000,F,[alpha0 alpha_free],'k--')
grid on
ylim([0 Inf])
xlabel('Bundle contraction (mm)')
ylabel('Bundle force (N)')
title('Design 3')
legend('Experiment','Model')

%% DESIGN 5
load('Mar13_13in_50psi_easier2plot.mat')
plot(DATA(:,2)*25.4,DATA(:,1)*1000*4.448)
hold on

P = 50*6894.76; % [Pa]
h0 = (1/16)/(3/8);
n = 1;
l0 = 13*0.0254;
beta0 = 0;
r0 = 0.3749*0.0254;
alpha0 = 30;

beta = @(alpha) real(asind((sind(beta0)*cosd(alpha0))/cosd(alpha)));
l = @(alpha) l0*(cosd(alpha)/cosd(alpha0));
C10 = 55.4298e3;
C01 = 260.5543e3;
t0 = h0*(r0);
B = l0/cosd(alpha0);     % [m] - length of fiber cord around FAM
N = (B*sind(alpha0))/(2*pi*r0); % [--] - number of fiber turns around FAM
F_gaylord = @(lambda) (P*((3*((lambda*l0)^2))-(B^2)))/(4*(N^2)*pi); 
Vb = (pi*(r0^2)*l0) - (pi*((r0-t0)^2)*l0); % [m^3] - bladder volume
LAMBDA = @(alpha) cosd(alpha)/cosd(alpha0);
term1 = @(lambda) (4*(C10+C01))*((l0^2)*(-1+(lambda^4)));
term2 = @(lambda) (4.*(l0.^6)*(-1+lambda)*(lambda^2)*(1+lambda)*(C10+(C01*(lambda^2))))/(((-4*(N^2)*(pi^2)*(r0^2))+((l0^2)*(-1+(lambda^2))))^2);
term3 = @(lambda) ((4*(l0^4)*(C10+(C01*(lambda^4))))/((-4*(N^2)*(pi^2)*(r0^2))+((l0^2)*(-1+(lambda^2)))));
term4 = @(lambda) ((l0^4)*(lambda^4)*(C10+(C01*(-1+(2*(lambda.^2))))))/((N^2)*(pi^2)*(r0^2));
F = @(lambda) F_gaylord(lambda) + Vb*((1./(2*(l0^3)*(lambda^3))))*(term1(lambda)+term2(lambda)-term3(lambda)-term4(lambda));
F = @(alpha) n*F(LAMBDA(alpha)).*cosd(beta(alpha)); % For each FAM pair in the bundle
alpha_range = alpha0:0.01:atand(sqrt(2));
for i = 1:length(alpha_range)
    model_deltalm2(i) = ((l0*cosd(beta0))-(l(alpha_range(i))*cosd(beta(alpha_range(i)))))/0.0254;
    model_F2(i) = F(alpha_range(i));
end
% alpha_c = fsolve(@(alpha) F(alpha),alpha0);
Fb = F(alpha0)*0.2248089431; % Blocked force (lbf)
alpha_free = fsolve(F,alpha0); % Braid angle at bundle free contraction (deg.)
deltalm_free = ((l0*cosd(beta0))-(l(alpha_free)*cosd(beta(alpha_free))))/0.0254; % Bundle free contraction 
% FAM rotation behavior changes due to inactive length
% l_in = (10.5*25.4)/1000;
% beta = @(alpha) asind((sind(beta0)*(l0+l_in))/(l(alpha)+l_in));
% fplot(@(alpha) (((l0+l_in)*cosd(beta0)) - ((l_in+l(alpha))*cosd(beta(alpha))))*1000,F,[alpha0 alpha_free],'k--')
% deltalm_free = (((l0+l_in)*cosd(beta0)) - ((l_in+l(alpha_free))*cosd(beta(alpha_free))))*1000; % Bundle free contraction [mm]
fplot(@(alpha) (((l0)*cosd(beta0)) - ((l(alpha))*cosd(beta(alpha))))*1000,F,[alpha0 alpha_free],'k--')
grid on
ylim([0 Inf])
xlabel('Bundle contraction (mm)')
ylabel('Bundle force (N)')
title('Design 5')
legend('Experiment','Model')

%% DESIGN 6
load('Mar8_14in_50psi_easier2plot.mat')
plot(DATA(:,2)*25.4,DATA(:,1)*1000*4.448)
hold on

P = 50*6894.76; % [Pa]
h0 = (1/16)/(3/8);
n = 1;
l0 = 14*0.0254;
beta0 = 0;
r0 = 0.3749*0.0254;
alpha0 = 30;

beta = @(alpha) real(asind((sind(beta0)*cosd(alpha0))/cosd(alpha)));
l = @(alpha) l0*(cosd(alpha)/cosd(alpha0));
C10 = 55.4298e3;
C01 = 260.5543e3;
t0 = h0*(r0);
B = l0/cosd(alpha0);     % [m] - length of fiber cord around FAM
N = (B*sind(alpha0))/(2*pi*r0); % [--] - number of fiber turns around FAM
F_gaylord = @(lambda) (P*((3*((lambda*l0)^2))-(B^2)))/(4*(N^2)*pi); 
Vb = (pi*(r0^2)*l0) - (pi*((r0-t0)^2)*l0); % [m^3] - bladder volume
LAMBDA = @(alpha) cosd(alpha)/cosd(alpha0);
term1 = @(lambda) (4*(C10+C01))*((l0^2)*(-1+(lambda^4)));
term2 = @(lambda) (4.*(l0.^6)*(-1+lambda)*(lambda^2)*(1+lambda)*(C10+(C01*(lambda^2))))/(((-4*(N^2)*(pi^2)*(r0^2))+((l0^2)*(-1+(lambda^2))))^2);
term3 = @(lambda) ((4*(l0^4)*(C10+(C01*(lambda^4))))/((-4*(N^2)*(pi^2)*(r0^2))+((l0^2)*(-1+(lambda^2)))));
term4 = @(lambda) ((l0^4)*(lambda^4)*(C10+(C01*(-1+(2*(lambda.^2))))))/((N^2)*(pi^2)*(r0^2));
F = @(lambda) F_gaylord(lambda) + Vb*((1./(2*(l0^3)*(lambda^3))))*(term1(lambda)+term2(lambda)-term3(lambda)-term4(lambda));
F = @(alpha) n*F(LAMBDA(alpha)).*cosd(beta(alpha)); % For each FAM pair in the bundle
alpha_range = alpha0:0.01:atand(sqrt(2));
for i = 1:length(alpha_range)
    model_deltalm2(i) = ((l0*cosd(beta0))-(l(alpha_range(i))*cosd(beta(alpha_range(i)))))/0.0254;
    model_F2(i) = F(alpha_range(i));
end
% alpha_c = fsolve(@(alpha) F(alpha),alpha0);
Fb = F(alpha0)*0.2248089431; % Blocked force (lbf)
alpha_free = fsolve(F,alpha0); % Braid angle at bundle free contraction (deg.)
deltalm_free = ((l0*cosd(beta0))-(l(alpha_free)*cosd(beta(alpha_free))))/0.0254; % Bundle free contraction 
% FAM rotation behavior changes due to inactive length
% l_in = (10.5*25.4)/1000;
% beta = @(alpha) asind((sind(beta0)*(l0+l_in))/(l(alpha)+l_in));
% fplot(@(alpha) (((l0+l_in)*cosd(beta0)) - ((l_in+l(alpha))*cosd(beta(alpha))))*1000,F,[alpha0 alpha_free],'k--')
% deltalm_free = (((l0+l_in)*cosd(beta0)) - ((l_in+l(alpha_free))*cosd(beta(alpha_free))))*1000; % Bundle free contraction [mm]
fplot(@(alpha) (((l0)*cosd(beta0)) - ((l(alpha))*cosd(beta(alpha))))*1000,F,[alpha0 alpha_free],'k--')
grid on
ylim([0 Inf])
xlabel('Bundle contraction (mm)')
ylabel('Bundle force (N)')
title('Design 6')
legend('Experiment','Model')
%% DESIGN 7
load('Mar13_16.7in_50psi_easier2plot.mat')
plot(DATA(:,2)*25.4,DATA(:,1)*1000*4.448)
hold on

P = 50*6894.76; % [Pa]
h0 = (1/16)/(3/8);
n = 2;
l0 = 16.75*0.0254;
beta0 = 6;
r0 = 0.3749*0.0254;
alpha0 = 30;

beta = @(alpha) real(asind((sind(beta0)*cosd(alpha0))/cosd(alpha)));
l = @(alpha) l0*(cosd(alpha)/cosd(alpha0));
C10 = 55.4298e3;
C01 = 260.5543e3;
t0 = h0*(r0);
B = l0/cosd(alpha0);     % [m] - length of fiber cord around FAM
N = (B*sind(alpha0))/(2*pi*r0); % [--] - number of fiber turns around FAM
F_gaylord = @(lambda) (P*((3*((lambda*l0)^2))-(B^2)))/(4*(N^2)*pi); 
Vb = (pi*(r0^2)*l0) - (pi*((r0-t0)^2)*l0); % [m^3] - bladder volume
LAMBDA = @(alpha) cosd(alpha)/cosd(alpha0);
term1 = @(lambda) (4*(C10+C01))*((l0^2)*(-1+(lambda^4)));
term2 = @(lambda) (4.*(l0.^6)*(-1+lambda)*(lambda^2)*(1+lambda)*(C10+(C01*(lambda^2))))/(((-4*(N^2)*(pi^2)*(r0^2))+((l0^2)*(-1+(lambda^2))))^2);
term3 = @(lambda) ((4*(l0^4)*(C10+(C01*(lambda^4))))/((-4*(N^2)*(pi^2)*(r0^2))+((l0^2)*(-1+(lambda^2)))));
term4 = @(lambda) ((l0^4)*(lambda^4)*(C10+(C01*(-1+(2*(lambda.^2))))))/((N^2)*(pi^2)*(r0^2));
F = @(lambda) F_gaylord(lambda) + Vb*((1./(2*(l0^3)*(lambda^3))))*(term1(lambda)+term2(lambda)-term3(lambda)-term4(lambda));
F = @(alpha) n*F(LAMBDA(alpha)).*cosd(beta(alpha)); % For each FAM pair in the bundle
alpha_range = alpha0:0.01:atand(sqrt(2));
for i = 1:length(alpha_range)
    model_deltalm2(i) = ((l0*cosd(beta0))-(l(alpha_range(i))*cosd(beta(alpha_range(i)))))/0.0254;
    model_F2(i) = F(alpha_range(i));
end
% alpha_c = fsolve(@(alpha) F(alpha),alpha0);
Fb = F(alpha0)*0.2248089431; % Blocked force (lbf)
alpha_free = fsolve(F,alpha0); % Braid angle at bundle free contraction (deg.)
% FAM rotation behavior changes due to inactive length
l_in = (10.5*25.4)/1000;
beta = @(alpha) asind((sind(beta0)*(l0+l_in))/(l(alpha)+l_in));
deltalm_free = (((l0+l_in)*cosd(beta0)) - ((l_in+l(alpha_free))*cosd(beta(alpha_free))))*1000; % Bundle free contraction [mm]
fplot(@(alpha) (((l0+l_in)*cosd(beta0)) - ((l_in+l(alpha))*cosd(beta(alpha))))*1000,F,[alpha0 alpha_free],'k--')
grid on
ylim([0 Inf])
xlabel('Bundle contraction (mm)')
ylabel('Bundle force (N)')
title('Design 7')
legend('Experiment','Model')