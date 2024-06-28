function mutation_offspring = mutation_homo_v2(xv,yv,L,W,D,alpha0,crossover_offspring)
    % INPUT ARGUMENTS: 
    % xv = vector of x-axis edge position of spatial envelope 
    % yv = vector of y-axis edge position of spatial envelope
    % L = length dimension of spatial envelope
    % W = width dimension of spatial envelope
    % D = depth dimension of spatial envelope 
    % alpha0 = initial braid angle
    % crossover_basepair = basepair properties after crossover [r0 beta0 l0
    % x_center y_center] 
    % pick_mutation = value indicating which property of base pair is selected
    % for mutation
    
    % OUTPUT: 
    % mutation_basepair = [r0 beta0 l0 x_center y_center]
    crossover_basepair = crossover_offspring(1,1:3);
    mutate_basepair = crossover_basepair;
%     num_mutate = randi(numel([1 2]));
    num_mutate = 1;
%     in = 1;
    if num_mutate == 1
%         while in == 1
            mutate_param = randi(numel([1 2 3]));
            % while mutate_param == 2 && crossover_basepair(2) == 0
            %     mutate_param = randi(numel([1 2 3]));
            % end
%             mutate_param = 1;
            if mutate_param == 1 % MUTATE RADIUS
                l = @(alpha) (cosd(alpha)/cosd(alpha0))*mutate_basepair(3);
                beta = @(alpha) real(asind((sind(mutate_basepair(2))*cosd(alpha0))/cosd(alpha)));
                beta_free = beta(atand(sqrt(2)));
                if beta_free < 90
                    alpha_free = atand(sqrt(2));
                else
                    alpha_free = real(acosd(sind(mutate_basepair(2))*cosd(alpha0)));
                end
                % ENSURES MUTATED RADIUS DIFFERS FROM CROSSOVER RADIUS
                while mutate_basepair(1) == crossover_basepair(1)
                    f = 0.7 + (1.3-0.7)*rand(1); % Controlled mutation
                    while f == 1
                        f = 0.7 + (1.3-0.7)*rand(1);
                    end
                    mutate_basepair(1) = f*crossover_basepair(1);
                    % CHECK MUTATED RADIUS SATISFIES WIDTH CONSTRAINT 
                    r = @(alpha) (sind(alpha)/sind(alpha0))*mutate_basepair(1);
                    width = @(alpha) (2*r(alpha)*cosd(beta(alpha)))+(l(alpha)*sind(beta(alpha)));
                    [~,max_width] = fminbnd(@(alpha) -width(alpha),alpha0,alpha_free);
                    while -max_width > (W/2) 
                        f = 0.7 + (1.3-0.7)*rand(1); % Controlled mutation
                        while f == 1
                            f = 0.7 + (1.3-0.7)*rand(1);
                        end
                        mutate_basepair(1) = f*crossover_basepair(1);
                        r = @(alpha) (sind(alpha)/sind(alpha0))*mutate_basepair(1);
                        width = @(alpha) (2*r(alpha)*cosd(beta(alpha)))+(l(alpha)*sind(beta(alpha)));
                        [~,max_width] = fminbnd(@(alpha) -width(alpha),alpha0,alpha_free);
                    end
                    % CHECK MUTATED RADIUS SATISFIES LENGTH CONSTRAINT
                    length = @(alpha) (2*r(alpha)*sind(beta(alpha)))+(l(alpha)*cosd(beta(alpha)));
                    [~,max_length] = fminbnd(@(alpha) -length(alpha),alpha0,alpha_free);
                    while -max_length > L
                        f = 0.7 + (1.3-0.7)*rand(1); % Controlled mutation
                        while f == 1
                            f = 0.7 + (1.3-0.7)*rand(1);
                        end
                        mutate_basepair(1) = f*crossover_basepair(1);
                        r = @(alpha) (sind(alpha)/sind(alpha0))*mutate_basepair(1);
                        length = @(alpha) (2*r(alpha)*sind(beta(alpha)))+(l(alpha)*cosd(beta(alpha)));
                        [~,max_length] = fminbnd(@(alpha) -length(alpha),alpha0,alpha_free);
                    end
                    % CHECK MUTATED RADIUS SATISFIES DEPTH CONSTRAINT
                    r0_max = (D/2)*(sind(alpha0)/sind(alpha_free));
                    while mutate_basepair(1) > r0_max || mutate_basepair(1) < 0
                        f = 0.7 + (1.3-0.7)*rand(1); % Controlled mutation
                        while f == 1
                            f = 0.7 + (1.3-0.7)*rand(1);
                        end
                        mutate_basepair(1) = f*crossover_basepair(1);
                    end
                end
%                 in = 0;
            end
            if mutate_param == 2 % MUTATE PENNATION ANGLE
                % ENSURES MUTATED PENNATION ANGLE DIFFERS FROM CROSSOVER
                % PENNATION ANGLE
                % ALSO, MUTATED PENNATION ANGLE MUST BE BTW 0 AND 90 DEG.
                alpha_free = real(asind((D/2)*(sind(alpha0)/mutate_basepair(1))));
                if alpha_free <= atand(sqrt(2))
                    beta0_min = real(asind(cosd(alpha_free)/cosd(alpha0)));
%                     beta0_max = 90;
                else
%                     beta0_max = 90;
                    beta0_min = 0;
                end
%                 crossover_basepair
%                 beta0_min
%                 beta0_max
                r = @(alpha) (sind(alpha)/sind(alpha0))*mutate_basepair(1);
                l = @(alpha) (cosd(alpha)/cosd(alpha0))*mutate_basepair(3);
                % USE WIDTH CONSTRAINT TO FIND THE LARGEST PENNATION ANGLE
                % syms beta0min
                width = @(beta0min) (2*r(alpha0)*cosd(beta0min))+(l(alpha0)*sind(beta0min))-(W/2);
                opts = optimset('Diagnostics','off', 'Display','off');
                beta0_MAX = fsolve(width,[0 90],opts);
                beta0_MAX = min(beta0_MAX);
                % beta0_MAX = vpasolve(width,beta0min,[0 90]);
                if isempty(beta0_MAX) == 1
                    width = @(beta0) (2*r(alpha0)*cosd(beta0))+(l(alpha0)*sind(beta0))-(W/2);
                    beta0_MAX = fsolve(width,30);
                end
%                 beta0_MAX
                if beta0_MAX < beta0_min
%                     in = 1;
                    mutation_offspring = [];
                    return
                end
                while mutate_basepair(2) == crossover_basepair(2) || mutate_basepair(2) > beta0_MAX || mutate_basepair(2) < beta0_min
                    f = rand(1);
                    while f == 1
                        f = rand(1);
                    end
                    mutate_basepair(2) = beta0_MAX*f;

%                     beta = @(alpha) real(asind((sind(mutate_basepair(2))*cosd(alpha0))/cosd(alpha)));
%                     beta_free = beta(atand(sqrt(2)));
%                     if beta_free < 90
%                         alpha_free = atand(sqrt(2));
%                     else
%                         alpha_free = real(acosd(sind(mutate_basepair(2))*cosd(alpha0)));
%                     end
%                     in = 0;
                end
            end

            if mutate_param == 3 % MUTATE LENGTH
                r = @(alpha) (sind(alpha)/sind(alpha0))*mutate_basepair(1);
                beta = @(alpha) real(asind((sind(mutate_basepair(2))*cosd(alpha0))/cosd(alpha)));
                beta_free = beta(atand(sqrt(2)));
                if beta_free < 90
                    alpha_free = atand(sqrt(2));
                else
                    alpha_free = real(acosd(sind(mutate_basepair(2))*cosd(alpha0)));
                end
                % ENSURES MUTATED LENGTH DIFFERS FROM CROSSOVER LENGTH & IS > 0
                while mutate_basepair(3) == crossover_basepair(3) || mutate_basepair(3) <= 0
                    f = 0.8 + (1.2-0.8)*rand(1); % Controlled mutation
                    while f == 1
                        f = 0.8 + (1.2-0.8)*rand(1);
                    end
                    mutate_basepair(3) = f*crossover_basepair(3);
                    l = @(alpha) (cosd(alpha)/cosd(alpha0))*mutate_basepair(3);
                    % CHECK MUTATED LENGTH SATISFIES WIDTH CONSTRAINT
                    width = @(alpha) (2*r(alpha)*cosd(beta(alpha)))+(l(alpha)*sind(beta(alpha)));
                    [~,max_width] = fminbnd(@(alpha) -width(alpha),alpha0,alpha_free);
                    while -max_width > (W/2) 
                        f = 0.7 + (1.3-0.7)*rand(1); % Controlled mutation
                        while f == 1
                            f = 0.7 + (1.3-0.7)*rand(1);
                        end
                        mutate_basepair(3) = f*crossover_basepair(3);
                        l = @(alpha) (cosd(alpha)/cosd(alpha0))*mutate_basepair(3);
                        width = @(alpha) (2*r(alpha)*cosd(beta(alpha)))+(l(alpha)*sind(beta(alpha)));
                        [~,max_width] = fminbnd(@(alpha) -width(alpha),alpha0,alpha_free);
                    end
                    % CHECK MUTATED LENGTH SATISFIES LENGTH CONSTRAINT
                    length = @(alpha) (2*r(alpha)*sind(beta(alpha)))+(l(alpha)*cosd(beta(alpha)));
                    [~,max_length] = fminbnd(@(alpha) -length(alpha),alpha0,alpha_free);
                    while -max_length > L
                        f = 0.7 + (1.3-0.7)*rand(1); % Controlled mutation
                        while f == 1
                            f = 0.7 + (1.3-0.7)*rand(1);
                        end
                        mutate_basepair(3) = f*crossover_basepair(3);
                        l = @(alpha) (cosd(alpha)/cosd(alpha0))*mutate_basepair(3);
                        length = @(alpha) (2*r(alpha)*sind(beta(alpha)))+(l(alpha)*cosd(beta(alpha)));
                        [~,max_length] = fminbnd(@(alpha) -length(alpha),alpha0,alpha_free);
                    end                
                end
%                 in = 0;
            end
%         end
    end
%     else
%         mutate_param(1) = randi(numel([1 2 3]));
%         mutate_param(2) = randi(numel([1 2 3]));
%         while mutate_param(2) == mutate_param(1)
%             mutate_param(2) = randi(numel([1 2 3]));
%         end
%     end
    r = @(alpha) (sind(alpha)/sind(alpha0))*mutate_basepair(1);
    l = @(alpha) (cosd(alpha)/cosd(alpha0))*mutate_basepair(3);
    beta = @(alpha) real(asind((sind(mutate_basepair(2))*cosd(alpha0))/cosd(alpha)));
    beta_free = beta(atand(sqrt(2)));
    if beta_free < 90
        alpha_free = atand(sqrt(2));
    else
        alpha_free = real(acosd(sind(mutate_basepair(2))*cosd(alpha0)));
    end
    % FIND CORRESPONDING xc,yc coordinate position
    x_min = @(alpha) (r(alpha)*cosd(beta(alpha)))+(l(alpha)*sind(beta(alpha)));
    [~,XC_MAX] = fminbnd(@(alpha) -x_min(alpha),alpha0,alpha_free);
    xc_max = (W/2)-(-XC_MAX);
    
    x1 = @(alpha) xc_max-(r(alpha)*cosd(beta(alpha)));
    [~,x1_min] = fminbnd(x1,alpha0,alpha_free);
    if x1_min < 0
        mutation_offspring = [];
        return
    end
    y_min = @(alpha) (r(alpha)*sind(beta(alpha)));
    [~,YC_MAX] = fminbnd(@(alpha) -y_min(alpha),alpha0,alpha_free);
    yc_max = L-(-YC_MAX);
    
    y3 = @(alpha) yc_max - (r(alpha)*sind(beta(alpha))) - (l(alpha)*cosd(beta(alpha)));
    [~,y3_min] = fminbnd(y3,alpha0,alpha_free);
    if y3_min < 0 
        % in = 0;
        mutation_offspring = [];
        return
    end

    % SOVLE FOR MAXIMUM NUMBER OF FAMS in bundle
    deltax = @(alpha) (2*r(alpha))/cosd(beta(alpha));
    [~,deltax_max] = fminbnd(@(alpha) -deltax(alpha),alpha0,alpha_free);
    deltax_max = -deltax_max;
    n = @(alpha) 2*((((W/2)-(2*r(alpha)*cosd(beta(alpha)))-(l(alpha)*sind(beta(alpha))))/deltax_max)+1);
    [~,nmax] = fminbnd(@(alpha) -n(alpha),alpha0,alpha_free);
    n_max = floor(-nmax);
    if rem(n_max,2) == 1
        n_max = n_max-1;
    end
    x1 = @(alpha) xc_max - ((n_max/2)-1)*deltax_max-(r(alpha)*cosd(beta(alpha)));
     [~,x1_min] = fminbnd(x1,alpha0,alpha_free);
    while x1_min < 0
        n_max = n_max-2;
        x1 = @(alpha) xc_max - ((n_max/2)-1)*deltax_max-(r(alpha)*cosd(beta(alpha)));
        [~,x1_min] = fminbnd(x1,alpha0,alpha_free);
        if n_max < 2
            mutation_offspring = [];
            return
        end
    end
    y3 = @(alpha) yc_max-(r(alpha)*sind(beta(alpha)))-(l(alpha)*cosd(beta(alpha)));
    [~,y3_min] = fminbnd(y3,alpha0,alpha_free);
    if y3_min < 0
        in = 0;
        mutation_offspring = [];
        return
    end
    % n = @(alpha) 2*(((L-(2*r(alpha)*sind(beta(alpha)))-(l(alpha)*cosd(beta(alpha))))/deltay_max)+1);
    % [~,nmax] = fminbnd(@(alpha) -n(alpha),alpha0,alpha_free);
    % n_max = floor(-nmax);
    % % ENSURE EVEN INTEGER
    % if rem(n_max,2) == 1
    %     n_max = n_max-1;
    % end
    % y4 = @(alpha) yc_max-((n_max/2)-1)*deltay_max-(r(alpha)*sind(beta(alpha)))-(l(alpha)*cosd(beta(alpha)));
    % [~,y4_min] = fminbnd(y4,alpha0,alpha_free);
    % while y4_min < 0 && n_max > 2
    %     n_max = n_max-2;
    %     y4 = @(alpha) yc_max-((n_max/2)-1)*deltay_max-(r(alpha)*sind(beta(alpha)))-(l(alpha)*cosd(beta(alpha)));
    %     [~,y4_min] = fminbnd(y4,alpha0,alpha_free);
    % end
    if n_max < 2
        mutation_offspring = [];
        return
    end
    if n_max == 2 
        r0 = mutate_basepair(1);    beta0 = mutate_basepair(2);     l0 = mutate_basepair(3);
        xc = xc_max;    yc = yc_max;
        mutation_offspring = [r0 beta0 l0 xc yc];
        return
    end
    if n_max > 2
        mutation_offspring = zeros(n_max/2,5);
        for i = 1:(n_max/2)
            mutation_offspring(i,1) = mutate_basepair(1);
            mutation_offspring(i,2) = mutate_basepair(2);
            mutation_offspring(i,3) = mutate_basepair(3);
            mutation_offspring(i,4) = xc_max-(i-1)*deltax_max;
            mutation_offspring(i,5) = yc_max;
        end
    end
    
end