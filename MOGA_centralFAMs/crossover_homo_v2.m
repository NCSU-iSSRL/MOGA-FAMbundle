% FUNCTION FOR DOING CROSSOVER (HOMOGENEOUS FAMs)
function crossover_offspring = crossover_homo_v2(xv,yv,L,W,D,alpha0,parents)
    % INPUT ARGUMENT:
    % parents = selected individuals in population
    
    % OUTPUTS:
    % crossover_offspring = properties after crossover in matrix form with n_max/2 rows for each FAM pair in the individual 
    
    num_parents = length(parents); % [--] - number of parents
    vector_parents = 1:num_parents; 
    
    % in = 0;
    % while in == 0
%         crossover_offspring = zeros(n_max/2,5);
        % RANDOMLY PICK TWO PARENTS TO CROSSOVER
        parent1 = randi(numel(vector_parents)); % [--] - select parent 1
        parent2 = randi(numel(vector_parents)); % [--] - select parent 2
        while parent1 == parent2 % Ensures parent 1 and parent 2 are NOT the same
            parent2 = randi(numel(vector_parents));
        end

        PARENT1 = parents{parent1}(1,:); % GENOME OF PARENT 1
        PARENT2 = parents{parent2}(1,:); % GENOME OF PARENT 2

        base_pair = zeros(1,3);

        % Choose which parent to get the pennation angle from 
        pick_beta0 = randi(numel([1 2]));
        if pick_beta0 == 1
            base_pair(2) = PARENT1(2);
        else
            base_pair(2) = PARENT2(2);
        end
        % Solve for corresponding pennation angle 
        beta = @(alpha) real(asind((sind(base_pair(2))*cosd(alpha0))/cosd(alpha)));
        beta_free = beta(atand(sqrt(2)));
        if base_pair(2) == 0
            beta = @(alpha) base_pair(2);
            beta_free = beta(atand(sqrt(2)));
        end
        if base_pair(2) == 90
            beta = @(alpha) base_pair(2);
            beta_free = base_pair(2);
        end
        if beta_free < 90
            alpha_free = atand(sqrt(2));
        else
            alpha_free = real(acosd(sind(base_pair(2))*cosd(alpha0)));
        end
        r0_max = (D/2)*(sind(alpha0)/sind(alpha_free));
        % Choose which parameter to pass down 
        options = [1 3];
        pick = randi(numel([1 2]));
        pass = options(pick);
        if pass == 1 % r0 is selected
            % Ensure r0 is selected from other parent
            if pick_beta0 == 1
                pick_r0 = 2;
            else
                pick_r0 = 1;
            end
            % Retrieve r0 from the parent
            if pick_r0 == 1
                base_pair(1) = PARENT1(1);
            else
                base_pair(1) = PARENT2(1);
            end
            r = @(alpha) (sind(alpha)/sind(alpha0))*base_pair(1);
            % Check selected r <= (D/2)
            if r(alpha_free) > r0_max
               base_pair(1) = r0_max;
            end
            % 2. CALCULATE MAXIMUM FAM LENGTH POSSIBLE WITH PENNATION ANGLE AND WIDTH DIMENSION OF SPATIAL ENVELOPE 
            l = @(alpha) (cosd(alpha0)/cosd(alpha))*(((W/2)-(2*r(alpha)*cosd(beta(alpha))))/sind(beta(alpha))); % [in] - function of FAM length 
            l0_max = l(fminbnd(l,alpha0,alpha_free)); % [in] - maximum possible FAM length 
            base_pair(3) = l0_max*rand(1); 
            while base_pair(3) == 0
                base_pair(3) = l0_max*rand(1); 
            end
            % Check FAM selected FAM length satisfies length dimension of spatial
            % envelope
            l_L = @(alpha) (2*r(alpha)*sind(beta(alpha)))+(cosd(alpha)/cosd(alpha0))*base_pair(3)*cosd(beta(alpha));
            [~,MAX] = fminbnd(@(alpha) -l_L(alpha),alpha0,alpha_free);
            while -MAX > L % If selected FAM length exceeds length dimension of spatial bounds, decrease length of FAM
                base_pair(3) = base_pair(3)*rand(1); % [in] - re-guess try another initial FAM length
                l_L = @(alpha) (2*r(alpha)*sind(beta(alpha)))+(cosd(alpha)/cosd(alpha0))*base_pair(3)*cosd(beta(alpha));
                [~,MAX] = fminbnd(@(alpha) -l_L(alpha),alpha0,alpha_free);
            end
%             l = @(alpha) (cosd(alpha)/cosd(alpha0))*base_pair(3);
        end
        if pass == 3 % l0 is selected
            % Ensure l0 is selected from other parent
            if pick_beta0 == 1
                pick_l0 = 2;
            else
                pick_l0 = 1;
            end
            % Retrieve l0 from the parent
            if pick_l0 == 1
                base_pair(3) = PARENT1(3);
            else
                base_pair(3) = PARENT2(3);
            end
            r = @(alpha) (sind(alpha)/sind(alpha0))*r0_max;
            % 2. CALCULATE MAXIMUM FAM LENGTH POSSIBLE WITH PENNATION ANGLE AND WIDTH DIMENSION OF SPATIAL ENVELOPE 
            l = @(alpha) (cosd(alpha0)/cosd(alpha))*(((W/2)-(2*r(alpha)*cosd(beta(alpha))))/sind(beta(alpha))); % [in] - function of FAM length 
            l0_max = l(fminbnd(l,alpha0,alpha_free)); % [in] - maximum possible FAM length 
            % Check selected l0 satisfies (2rcos(beta))+(lsin(beta))<=(W/2)
            if base_pair(3) > l0_max
                base_pair(3) = l0_max;
            end
            % Check FAM selected FAM length satisfies length dimension of spatial
            % envelope
            l_L = @(alpha) (2*r(alpha)*sind(beta(alpha)))+(cosd(alpha)/cosd(alpha0))*base_pair(3)*cosd(beta(alpha));
            [~,MAX] = fminbnd(@(alpha) -l_L(alpha),alpha0,alpha_free);
            if -MAX > L % If selected FAM length exceeds length dimension of spatial bounds, decrease length of FAM
                l_L = @(alpha) ((L-(2*r(alpha)*sind(beta(alpha))))/cosd(beta(alpha)));
                [~,MAX] = fminbnd(@(alpha) -l_L(alpha),alpha0,alpha_free);
                base_pair(3) = -MAX; 
            end
            l = @(alpha) (cosd(alpha)/cosd(alpha0))*base_pair(3);
            r_W = @(alpha) ((W/2)-(l(alpha)*sind(beta(alpha))))/(2*cosd(beta(alpha)));
            [~,R0_MAX] = fminbnd(r_W,alpha0,alpha_free);
            base_pair(1) = 0.1 + R0_MAX*rand(1);
            while base_pair(1) == 0
                base_pair(1) = 0.1 + R0_MAX*rand(1); 
            end
            if base_pair(1) > r0_max
                base_pair(1) = r0_max;
            end
        end
        r = @(alpha) (sind(alpha)/sind(alpha0))*base_pair(1);
        l = @(alpha) (cosd(alpha)/cosd(alpha0))*base_pair(3);

        x_min = @(alpha) (r(alpha)*cosd(beta(alpha)))+(l(alpha)*sind(beta(alpha)));
        [~,XC_MAX] = fminbnd(@(alpha) -x_min(alpha),alpha0,alpha_free);
        xc_max = (W/2)-(-XC_MAX);
        
        x1 = @(alpha) xc_max-(r(alpha)*cosd(beta(alpha)));
        [~,x1_min] = fminbnd(x1,alpha0,alpha_free);
        if x1_min < 0
            % in = 0;
            crossover_offspring = [];
            return
        end
        
        y_min = @(alpha) (r(alpha)*sind(beta(alpha)));
        [~,YC_MAX] = fminbnd(@(alpha) -y_min(alpha),alpha0,alpha_free);
        yc_max = L-(-YC_MAX);

        deltay = @(alpha) (2*r(alpha))/sind(beta(alpha));
        [~,deltay_max] = fminbnd(@(alpha) -deltay(alpha),alpha0,alpha_free);
        deltay_max = -deltay_max;
        n = @(alpha) 2*(((L-(2*r(alpha)*sind(beta(alpha)))-(l(alpha)*cosd(beta(alpha))))/deltay_max)+1);
        [~,nmax] = fminbnd(@(alpha) -n(alpha),alpha0,alpha_free);
        n_max = floor(-nmax);
        if rem(n_max,2) == 1
            n_max = n_max-1;
        end
        y4 = @(alpha) yc_max-((n_max/2)-1)*deltay_max-(r(alpha)*sind(beta(alpha)))-(l(alpha)*cosd(beta(alpha)));
        [~,y4_min] = fminbnd(y4,alpha0,alpha_free);
        while y4_min < 0
            n_max = n_max-2;
            y4 = @(alpha) yc_max-((n_max/2)-1)*deltay_max-(r(alpha)*sind(beta(alpha)))-(l(alpha)*cosd(beta(alpha)));
            [~,y4_min] = fminbnd(y4,alpha0,alpha_free);
        end
        % n = @(alpha) 2*(((L-(2*r(alpha)*sind(beta(alpha)))-(l(alpha)*cosd(beta(alpha))))/deltay_max)+1);
        % [~,nmax] = fminbnd(@(alpha) -n(alpha),alpha0,alpha_free);
        % n_max = floor(-nmax);
        % if rem(n_max,2) == 1
        %     n_max = n_max-1;
        % end
        % y4 = @(alpha) yc_max-((n_max/2)-1)*deltay_max-(r(alpha)*sind(beta(alpha)))-(l(alpha)*cosd(beta(alpha)));
        % [~,y4_min] = fminbnd(y4,alpha0,alpha_free);
        % while y4_min < 0 && in == 1 && n_max >= 2
        %     n_max = n_max-2;
        %     y4 = @(alpha) yc_max-((n_max/2)-1)*deltay_max-(r(alpha)*sind(beta(alpha)))-(l(alpha)*cosd(beta(alpha)));
        %     [~,y4_min] = fminbnd(y4,alpha0,alpha_free);
        % end
        if isnan(n_max) == 1 || isnan(base_pair(1)) == 1
%             in = 0;
            crossover_offspring = [];
            return
        end
        if y4_min < 0 && n_max < 2
%             in = 0;
            crossover_offspring = [];
            return
        else 
            % in = 1;
            crossover_offspring = zeros(n_max/2,5);
            for i = 1:n_max/2
                crossover_offspring(i,1) = base_pair(1);
                crossover_offspring(i,2) = base_pair(2);
                crossover_offspring(i,3) = base_pair(3);
                crossover_offspring(i,4) = xc_max;
                crossover_offspring(i,5) = yc_max-(i-1)*deltay_max;
            end
            return
        end
    % end

%     n_max
%     base_pair
%     crossover_offspring = zeros(n_max/2,5);

end
