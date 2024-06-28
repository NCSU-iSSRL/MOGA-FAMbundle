clc
clear all
close all

% Establish spatial envelope of interest
L = 24;     W = 12;  D = 2;    % Length, width, and depth dimension of spatial envelope
xv = [0 W W 0 0];   yv = [0 0 L L 0]; % For a rectangle
envelope = polyshape(xv,yv);
[xv,yv] = boundary(envelope);
% xv = [0 W 0.75*W 0.25*W 0]; yv = [0 0 L L 0];   % For a trapezoid
% xv = [0 W 0.75*W 0.25*W 0];   yv = [L L 0 0 L]; % For an inverted trapezoid
% For an ellipse
% a = (W/2);  b = (L/2);  x0 = a;     y0 = b;     t = -pi:0.01:pi;
% xv = x0+a*cos(t);   yv = y0+b*sin(t);
alpha0 = 30;

% GENERATE INDIVIDUALS FOR INITAL POPULATION (100 individuals) --> (200 individuals)
num_individuals = 200;
population = cell(1,num_individuals);
figure
tic
for I = 1:num_individuals
    I
    population{1,I} = generate_individual_homo_v2(xv,yv,L,W,D,alpha0);
    while isempty(population{1,I}) == 1
        population{1,I} = generate_individual_homo_v2(xv,yv,L,W,D,alpha0);
    end
    toc
    pause(0.01)
end
save('homopop_200_rect_central_aug16.mat')

%%
clc
clear all
close all
% load('maxmusclevolume_pinned_Fb_N.mat')
% load('maxmusclevolume_pinned_Vm_cm3.mat')
load('homopop_200_rect_central_aug16.mat')
% for i = 1:num_individuals
%     [overlap(i),~] = check_overlap_obj(xv,yv,L,W,D,alpha0,population{i});
%     plot(i,overlap(i),'k.')
%     hold on
%     pause(0.01)
% end
fitness = calc_fitness_Fb_deltalm(xv,yv,L,W,D,alpha0,population);
save('homopop_200_rect_central_aug16_Fb_deltalm.mat','fitness')