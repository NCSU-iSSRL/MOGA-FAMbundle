% FUNCTION TO CALCULATE FITNESS OF EACH INDIVIDUAL IN POPULATION
function fitness = calc_fitness_Fb_kh_deltalm(xv,yv,L,W,D,alpha0,population)
% function [fitness1, fitness2] = calc_fitness(xv,yv,L,W,D,alpha0,population)
    % INPUT ARGUMENTS:
    % xv = vector of x-axis edge position of spatial envelope 
    % yv = vector of y-axis edge position of spatial envelope
    % L = length dimension of spatial envelope
    % W = width dimension of spatial envelope
    % D = depth dimension of spatial envelope 
    % alpha0 = initial braid angle
    % population = cell array of individuals in population
    % NOTE: r0 [in] beta0 [deg] l0 [in] x_center [in] y_center [in]
    
    % OUTPUT ARGUMENTS: 
    % fitness = fitness value of each individual in population 

    num_individuals = length(population);
    % subplot(2,2,4)
%     C10 = 62.9e3;       C01 = 113.9e3; % [Pa] - Mooney Rivlin material constants for silicone bladder 
%     h0 = 0.1; % [--] - bladder wall thickness ratio (t0/r0)
%     P = 100*6894.76; % [Pa] - applied pressure to FAM
%     LAMBDA = @(alpha) (cosd(alpha)/cosd(alpha0));
%     deltay = zeros(1,num_individuals);
    fitness = zeros(num_individuals,2);
    for i = 1:num_individuals
%         i
        F_total = @(alpha) 0;
        r0 = zeros(1,length(population{i}(:,1)));
        beta0 = zeros(1,length(population{i}(:,1)));
        l0 = zeros(1,length(population{i}(:,1)));
%         x_center = zeros(1,length(population{i}(:,1)));
        y_center = zeros(1,length(population{i}(:,1)));
        % r = cell(1,length(population{i}(:,1)));
        % beta = cell(1,length(population{i}(:,1)));
        % l = cell(1,length(population{i}(:,1)));
        % h0 = 0.1;
        h0 = (1/16)/(3/8);
        t0 = zeros(1,length(population{i}(:,1)));
%         alpha_free = zeros(1,length(population{i}(:,1)));
        F = cell(1,length(population{i}(:,1)));
%         V = zeros(1,length(population{i}(:,1)));
        % i
        for j = 1:length(population{i}(:,1))
            % population{6}(j,1)*0.0254
            tf = isa(population{i}(j,:),'cell');
            if tf == 1
                FAM = cell2mat(population{i}(j,:));
            else
                FAM = cell2mat(num2cell(population{i}(j,:)));
            end
            % FAM = population{i};
            % Convert dimensions all from inches to meters
            r0(j) = FAM(1)*0.0254;      beta0(j) = FAM(2);   l0(j) = FAM(3)*0.0254;
%             x_center(j) = FAM(4)*0.0254;    
            y_center(j) = FAM(5)*0.0254;
            % r0(j) = population{i}{1,1}(j,1)*0.0254;    beta0(j) = population{i}{1,1}(j,2);     l0(j) = population{i}{1,1}(j,3)*0.0254;
            % x_center(j) = population{i}{1,1}(j,4).*0.0254;  
            y_center(j) = population{i}(j,5).*0.0254;

%             V(j) = (2*pi*(r0(j)^2)*l0(j));

            % r = @(alpha) (sind(alpha)./sind(alpha0))*r0(j);
            % l = @(alpha) (cosd(alpha)./cosd(alpha0))*l0(j);
            beta = @(alpha) real(asind((sind(beta0(j))*cosd(alpha0))./cosd(alpha)));
%             y{j} = @(alpha) y_center(j)-l(alpha).*cosd(beta(alpha));
        
            t0(j) = h0*(r0(j));
            P = 50*6894.76;
            % C10 = 62.9e3;
            % C01 = 113.9e3;
            C10 = 55.4298e3;
            % C01 = 550.5543e3;
            C01 = 260.5543e3;
            B = l0(j)/cosd(alpha0);     % [m] - length of fiber cord around FAM
            N = (B*sind(alpha0))/(2*pi*r0(j)); % [--] - number of fiber turns around FAM
            F_gaylord = @(lambda) (P*((3*((lambda*l0(j)).^2))-(B^2)))/(4*(N^2)*pi); 
            Vb = (pi*(r0(j)^2)*l0(j)) - (pi*((r0(j)-t0(j))^2)*l0(j)); % [m^3] - bladder volume
            LAMBDA = @(alpha) cosd(alpha)/cosd(alpha0);
            term1 = @(lambda) (4*(C10+C01)).*((l0(j)^2).*(-1+(lambda.^4)));
            term2 = @(lambda) (4.*(l0(j).^6).*(-1+lambda).*(lambda.^2).*(1+lambda).*(C10+(C01.*(lambda.^2))))./(((-4*(N^2)*(pi^2)*(r0(j)^2))+((l0(j)^2).*(-1+(lambda.^2)))).^2);
            term3 = @(lambda) ((4*(l0(j)^4)*(C10+(C01.*(lambda.^4))))./((-4*(N^2)*(pi^2)*(r0(j)^2))+((l0(j)^2).*(-1+(lambda.^2)))));
            term4 = @(lambda) ((l0(j)^4).*(lambda.^4).*(C10+(C01.*(-1+(2.*(lambda.^2))))))/((N^2)*(pi^2)*(r0(j)^2));
            F{j} = @(lambda) F_gaylord(lambda) + Vb.*((1./(2*(l0(j)^3).*(lambda.^3)))).*(term1(lambda)+term2(lambda)-term3(lambda)-term4(lambda));
            F{j} = @(alpha) 2*F{j}(LAMBDA(alpha)).*cosd(beta(alpha)); % For each FAM pair in the bundle
            % F{j} = @(alpha) 2*F_gaylord(LAMBDA(alpha)).*cosd(beta(alpha));
            % fplot(F{j})
            % hold on
            % pause(0.1)
            % F{j}(alpha0)
            % fzero(F{j},alpha0)
            % alpha_free(j) = fzero(F{j},alpha0);
            F_total = @(alpha) F_total(alpha) + F{j}(alpha);
%             fplot(y{j},[alpha0 alpha_free(j)])
%             fplot(F{j},[alpha0 alpha_free(j)])
%             hold on
        end
        F_block = F_total(alpha0);
        beta_free = beta(atand(sqrt(2)));
        % i
        if beta_free < 90
            alpha_free = atand(sqrt(2));
        else
            alpha_free = real(acosd(sind(beta0(end))*cosd(alpha0)));
        end
        % alpha_free
        lmi = y_center(end)-(l0(end)*cosd(beta0(end)));
        r = @(alpha) (sind(alpha)./sind(alpha0))*r0(end);
        l = @(alpha) (cosd(alpha)./cosd(alpha0))*l0(end);
        lmf = y_center(end)-(l(alpha_free)*cosd(beta(alpha_free)));
        % deltay = lmf-lmi;
        options = optimset('display','off');
        [alpha_fzero,~] = fsolve(F_total,alpha0,options);
        deltay = @(alpha) (y_center(end)-(l(alpha)*cosd(beta(alpha))))-lmi;
        deltalm_free = deltay(alpha_fzero);
        % alpha_range = linspace(alpha0,alpha_free,10e3);
        % F_TOTAL = F_total(alpha_range);
        % Index = find(F_TOTAL<0);
        % if isempty(Index) == 1
        %     deltay = lmf-lmi;
        % else
        %     delta = @(alpha) lmi-(y_center(end)-(l(alpha).*cosd(beta(alpha))));
        %     deltay = delta(alpha_range(Index(1)-1));
        % end
        % if deltay < 0 
        %     F_block = 0;
        %     deltay = 0;
        % end
%         tf = isa(population{i},'cell');
%         if tf == 0
%             population{i} = population(i);
%             % population{i} = mat2cell(population{i},[1]);
%         end
% %         population{i}
% %         population{i}{1,1}
%         [overlap,alpha_max] = check_overlap_obj(xv,yv,L,W,D,alpha0,population{i}{1,1});
%         if overlap == 1
%             deltay(i) = 0;
% %             F_block = 0;
%         end
%         if overlap == 0
%             % alpha_interval = alpha0:0.1:alpha_max;
%             alpha_max = alpha_max(end);
%             l_max = l0(end)*(cosd(alpha_max)/cosd(alpha0));
%             beta_max = real(asind((sind(beta0(end))*cosd(alpha0))/cosd(alpha_max)));
%             deltay(i) = (y_center(end)-l_max*cosd(beta_max))-(y_center(end)-l0(end)*cosd(beta0(end)));
%             F_block = F_total(alpha0);
%         end
        fitness(i,1) = F_block; % [N] - blocked force
        fitness(i,2) = deltalm_free; % [m] - displacement
%         fitness(i,2) = sum(V);
        clear F_total r0 beta0 l0 x_center y_center r l t0 y F V F_block lmi lmf deltalm_free
    end
    
end