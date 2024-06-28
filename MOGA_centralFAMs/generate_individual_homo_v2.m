% FUNCTION FOR GENERATING A FEASIBLE BASE PAIR FOR CENTRALLY ARRANGED
% PENNATE FAM BUNDLE
function individual = generate_individual_homo_v2(xv,yv,L,W,D,alpha0)
    % xv = x position vector of spatial envelope
    % yv = y position vector of spatial envelope
    % L = length dimension of spatial envelope
    % W = width dimension of spatial envelope
    % D = depth dimension of spatial envelope
    % alpha0 = initial FAM braid angle
    % figure
    % plot(xv,yv,'k','LineWidth',1.5) % Display spatial envelope
    % xline(W/2,'k--') % Display axial centerline
    % hold on
    % axis equal 
    % xlim([0 W])
    % ylim([0 L])

    beta0_min = 0;  % [deg] - minimum possible FAM pennation angle
    beta0_max = 90; % [deg] - maximum possible FAM pennation angle
    
    % 1. RANDOMLY SELECT AN INITIAL PENNATION ANGLE FOR THE FAM 
    beta0 = beta0_min + (beta0_max-beta0_min)*rand(1); % [deg] - randomly picking an initial pennation angle 
    beta = @(alpha) real(asind((sind(beta0)*cosd(alpha0))./cosd(alpha))); % [deg] - function for instantaneous FAM pennation angle 
    beta_free = beta(atand(sqrt(2))); % [deg] - FAM pennation angle at free contraction
    if beta0 == 0
        beta = @(alpha) beta0;
        beta_free = beta(atand(sqrt(2)));
    end
    if beta0 == 90
        r0 = d0/2;
        l0 = 3;
        beta = @(alpha) beta0;
        beta_free = beta(atand(sqrt(2)));
        x_center = 0;
        y_center = L-r0;
        feasible = 1;
        base_pair = [r0 beta0 l0 x_center y_center];
        return
    end
    % Based on FAM pennation angle at free strain, calculate FAM braid angle at
    % free strain
    if beta_free == 90
        alpha_free = real(acosd(sind(beta0)*cosd(alpha0))); % [deg] - ROTATION-LIMITED
    else
        alpha_free = atand(sqrt(2)); % [deg] - CONTRACTION LIMITED
    end
    
    d0_min = D/2; % [in] - minimum possible initial outer diameter of FAM
    d0_max = D*(sind(alpha0)/sind(alpha_free)); % [in] - maximum possible initial outer diameter of FAM 
    d0 = d0_min + (d0_max-d0_min)*rand(1); % [in] - randomly picking an initial outer diameter of FAM
    r0 = d0/2; % [in] - initial outer radius of FAM
    r = @(alpha) (sind(alpha)./sind(alpha0))*r0; % [in] - function for instantaneous FAM outer radius
    
    l0_min = 0.1;
    % 2. CALCULATE MAXIMUM FAM LENGTH POSSIBLE WITH PENNATION ANGLE AND WIDTH DIMENSION OF SPATIAL ENVELOPE 
    l = @(alpha) (cosd(alpha0)/cosd(alpha))*(((W/2)-(2*r(alpha)*cosd(beta(alpha))))/sind(beta(alpha))); % [in] - function of FAM length 
    l0_max = l(fminbnd(l,alpha0,alpha_free)); % [in] - maximum possible FAM length 
    % 3. RANDOMLY SELECT INITIAL FAM LENGTH <= MAXIMUM FAM LENGTH POSSIBLE
    l0 = (l0_min) + (l0_max-l0_min)*rand(1); % [in] - randomly picking an initial FAM length
    % Check FAM selected FAM length satisfies length dimension of spatial
    % envelope
    l_L = @(alpha) (2*r(alpha)*sind(beta(alpha)))+(cosd(alpha)/cosd(alpha0))*l0*cosd(beta(alpha));
    [~,MAX] = fminbnd(@(alpha) -l_L(alpha),alpha0,alpha_free);
    while -MAX > L % If selected FAM length exceeds length dimension of spatial bounds, decrease length of FAM
        l0 = (l0_min) + (l0-l0_min)*rand(1); % [in] - re-guess try another initial FAM length
        l_L = @(alpha) (2*r(alpha)*sind(beta(alpha)))+(cosd(alpha)/cosd(alpha0))*l0*cosd(beta(alpha));
        [~,MAX] = fminbnd(@(alpha) -l_L(alpha),alpha0,alpha_free);
    end
    l = @(alpha) (cosd(alpha)/cosd(alpha0))*l0;

    % XC = @(alpha) (W/2)-(r(alpha)*cosd(beta(alpha)))-(l(alpha)*sind(beta(alpha)));
    % [~,xc_max] = fminbnd(@(alpha) XC(alpha),alpha0,alpha_free);
    % xc_max = xc_max;
    % XC = @(alpha) r(alpha)*cosd(beta(alpha));
    % [~,xc_min] = fminbnd(@(alpha) -XC(alpha),alpha0,alpha_free);
    % xc_min = -xc_min;
    % 
    % YC = @(alpha) L-(r(alpha)*sind(beta(alpha)));
    % [~,yc_max] = fminbnd(@(alpha) -YC(alpha),alpha0,alpha_free);
    % yc_max = -yc_max;
%     [in_c,on] = inpolygon(xc_max,yc_max,xv,yv);
    x_min = @(alpha) (r(alpha)*cosd(beta(alpha)))+(l(alpha)*sind(beta(alpha)));
    [~,XC_MAX] = fminbnd(@(alpha) -x_min(alpha),alpha0,alpha_free);
    xc_max = (W/2)-(-XC_MAX);
    
    y_min = @(alpha) (r(alpha)*sind(beta(alpha)));
    [~,YC_MAX] = fminbnd(@(alpha) -y_min(alpha),alpha0,alpha_free);
    yc_max = L-(-YC_MAX);

    % x1 = xc_max-r0*cosd(beta0);   y1 = yc_max-r0*sind(beta0);
    % x2 = xc_max+r0*cosd(beta0);   y2 = yc_max+r0*sind(beta0);
    % x3 = x2+l0*sind(beta0);         y3 = y2-l0*cosd(beta0);
    % x4 = x1+l0*sind(beta0);         y4 = y1-l0*cosd(beta0);
%     [in,on] = inpolygon([x1 x2 x3 x4],[y1 y2 y3 y4],xv,yv);
%     individual = [r0 beta0 l0 xc_max yc_max];
%     [overlap,alpha_max] = check_overlap_obj(xv,yv,L,W,D,alpha0,individual);
    % patch([x1 x2 x3 x4],[y1 y2 y3 y4],'r','FaceAlpha',.3)
    % pause(0.1)

    deltay = @(alpha) (2*r(alpha))/sind(beta(alpha));
    [~,deltay_max] = fminbnd(@(alpha) -deltay(alpha),alpha0,alpha_free);
    deltay_max = -deltay_max;
    % n_max = 2*(((L-yc_max-(r0*sind(beta0))-(l0*cosd(beta0)))*deltay_max)+1)
    n = @(alpha) 2*(((L-(2*r(alpha)*sind(beta(alpha)))-(l(alpha)*cosd(beta(alpha))))/deltay_max)+1);
    [~,nmax] = fminbnd(@(alpha) -n(alpha),alpha0,alpha_free);
    n_max = floor(-nmax);
    if rem(n_max,2) == 1
        n_max = n_max-1;
    end
            
    % for i = 1:(n_max/2)
    %     x1 = xc_max-r0*cosd(beta0);   y1 = yc_max-(i-1)*deltay_max-r0*sind(beta0);
    %     x2 = xc_max+r0*cosd(beta0);   y2 = yc_max-(i-1)*deltay_max+r0*sind(beta0);
    %     x3 = x2+l0*sind(beta0);         y3 = y2-l0*cosd(beta0);
    %     x4 = x1+l0*sind(beta0);         y4 = y1-l0*cosd(beta0);
    %     patch([x1 x2 x3 x4],[y1 y2 y3 y4],'b','FaceAlpha',.3)
    %     hold on
    % end
    y4 = @(alpha) yc_max-((n_max/2)-1)*deltay_max-(r(alpha)*sind(beta(alpha)))-(l(alpha)*cosd(beta(alpha)));
    [~,y4_min] = fminbnd(y4,alpha0,alpha_free);
    while y4_min < 0
        n_max = n_max-2;
        y4 = @(alpha) yc_max-((n_max/2)-1)*deltay_max-(r(alpha)*sind(beta(alpha)))-(l(alpha)*cosd(beta(alpha)));
        [~,y4_min] = fminbnd(y4,alpha0,alpha_free);
    end
    % y4 = yc_max - ((n_max/2)-1)*deltay_max-(r(alpha0)*sind(beta(alpha0)))-(l(alpha0)*cosd(beta(alpha0)));
    % while y4 < 0
    %     n_max = n_max-2;
    %     y4 = yc_max - ((n_max/2)-1)*deltay_max-(r(alpha0)*sind(beta(alpha0)))-(l(alpha0)*cosd(beta(alpha0)));
    % end
    for i = 1:(n_max/2)
        x1 = xc_max-r0*cosd(beta0);   y1 = yc_max-(i-1)*deltay_max-r0*sind(beta0);
        x2 = xc_max+r0*cosd(beta0);   y2 = yc_max-(i-1)*deltay_max+r0*sind(beta0);
        x3 = x2+l0*sind(beta0);         y3 = y2-l0*cosd(beta0);
        x4 = x1+l0*sind(beta0);         y4 = y1-l0*cosd(beta0);
        % patch([x1 x2 x3 x4],[y1 y2 y3 y4],'r','FaceAlpha',.3)
        % hold on
        % patch(W-[x1 x2 x3 x4],[y1 y2 y3 y4],'r','FaceAlpha',.3)
        % hold on
    end
    individual = zeros(n_max/2,5);
    for i = 1:n_max/2
        individual(i,1) = r0;
        individual(i,2) = beta0;
        individual(i,3) = l0;
        individual(i,4) = xc_max;
        individual(i,5) = yc_max-(i-1)*deltay_max;
    end
%     [overlap,~] = check_overlap_obj(xv,yv,L,W,D,alpha0,individual)
end