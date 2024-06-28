% Emily Duan
clc
clear all
close all
% % SCRIPT FOR MULTIOBJECTIVE OPTIMIZATION GENETIC ALGORITHM % %
warning('off','all')
warning

% Establish spatial envelope of interest
D = 1;
L = 6*D;     W = 12*D;      % Length, width, and depth dimension of spatial envelope
xv = [0 W W 0 0];   yv = [0 0 L L 0]; % For a rectangle
envelope = polyshape(xv,yv);
[xv,yv] = boundary(envelope);
alpha0 = 30;

% GENERATE INDIVIDUALS FOR INITAL POPULATION (100 individuals) --> (200 individuals)
num_individuals = 200;
population = cell(1,num_individuals);
for I = 1:num_individuals
    population{1,I} = generate_individual_homo_v2(xv,yv,L,W,D,alpha0);
    while isempty(population{1,I}) == 1
        population{1,I} = generate_individual_homo_v2(xv,yv,L,W,D,alpha0);
    end
end
fitness = calc_fitness_Fb_kh_deltalm(xv,yv,L,W,D,alpha0,population);

num_individuals = length(population); % [--] - number of individuals in population 
num_gen = 5000; % [--] - number of generations 
num_objectives = length(fitness(1,:)); % [--] - number of objectives in optimization problem

% Initialize variables 
FITNESS = cell(1,num_gen); % fitness cell array for individuals in population at each generation
POPULATION = cell(num_gen,num_individuals); % population cell array for individuals in population at each generation

FITNESS{1} = fitness; 
% 3. NON-DOMINATING SORT FITNESS FROM INITIAL POPULATION 
sorted_I = ndsort(FITNESS{1},num_individuals); 
% Initialize variables
sorted_pop = cell(1,num_individuals);
sorted_fitness = zeros(num_individuals,num_objectives);
% Re-arrange individuals and their corresponding fitness in population based on sorting method 
for k = 1:num_individuals
    sorted_pop{k} = population{sorted_I(k)};
    sorted_fitness(k,:) = fitness(sorted_I(k),:);
end
POPULATION(1,:) = sorted_pop;
FITNESS{1} = sorted_fitness;
tic
for j = 1:num_gen
    j
    % SELECT ELITE
    num_elite = 0.3*num_individuals;
    elite = POPULATION(j,1:num_elite);
    % SELECT PARENTS
    num_parents = 0.5*num_individuals;
    parents = POPULATION(j,1:num_parents);
%     num_offspring = 1;
    num_mutate = 0.3*num_individuals;
    num_offspring = num_individuals-num_elite;
    offspring = cell(1,num_offspring);
    count = 1;
    for i = 1:num_offspring
        if i <= num_mutate
            offspring{i} = mutation_homo_v2(xv,yv,L,W,D,alpha0,parents{i}(1,:));
            while isempty(offspring{i}) == 1
                offspring{i} = mutation_homo_v2(xv,yv,L,W,D,alpha0,parents{i}(1,:));
            end
        else
            offspring{i} = crossover_homo_v2(xv,yv,L,W,D,alpha0,parents);
            while isempty(offspring{i}) == 1
                offspring{i} = crossover_homo_v2(xv,yv,L,W,D,alpha0,parents);
            end
        end
    end
    POPULATION(j+1,1:num_elite) = elite;
    POPULATION(j+1,num_elite+1:end) = offspring;
    for k = 1:num_individuals 
        POPULATION{j+1,k} = cell2mat(POPULATION(j+1,k));
    end
    FITNESS{j+1} = calc_fitness_Fb_kh_deltalm(xv,yv,L,W,D,alpha0,POPULATION(j+1,:));
    sorted_I = ndsort(FITNESS{j+1},num_individuals);
    sorted_pop = cell(1,num_individuals);
    for k = 1:num_individuals
        sorted_pop{k} = POPULATION{j+1,sorted_I(k)};
        sorted_fitness(k,:) = FITNESS{j+1}(sorted_I(k),:);
    end
    POPULATION(j+1,:) = sorted_pop;
    FITNESS{j+1} = sorted_fitness;
end
toc
save('maxFbkh_maxdeltalm_5000gen_rect_lateral_mar19_6x12x1_0.3E_0.5P_0.3M_0.1667h0_trial4_C01_260.5543e3.mat')