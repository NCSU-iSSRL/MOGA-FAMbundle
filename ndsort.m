% SCRIPT FOR NON-DOMINATED SORTING
function nondom_sort = ndsort(fitness,num_individuals)
    % INPUT: 
    % fitness = matrix of fitness values for each individual in population
    % num_individuals = number of individuals in population

    for i = 1:num_individuals
        dom_count(i) = 0;
        check_i = 1:num_individuals;
        check_i = check_i(check_i ~= i);
        for j = 1:length(check_i)
            % Depending on objectives: 
            % In this case, we are maximizing fitness(:,1) and maximizing
            % fitness(:,2)
            if (fitness(i,1) <= fitness(check_i(j),1) && fitness(i,2) <= fitness(check_i(j),2)) && (fitness(i,1) < fitness(check_i(j),1) || fitness(i,2) < fitness(check_i(j),2))
                dom_count(i) = dom_count(i)+1;
            end
        end
    end
    [~,I] = sort(dom_count,'ascend');
    % min(dom_count)
    % max(dom_count)
    nondom_sort = I;
    % nondom_sort = cell(1,num_individuals);
    % for i = 1:length(I)
    %     nondom_sort{i} = population{I(i)};
    % end
end