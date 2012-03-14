% Simplex Stochastic Collocation as described in
% "Refinement Criteria" SIAM Paper

function main()
    clear all; clc; 

    global TEST_FUNC;
    global p;

%-------------------------------------------------------------------------%
%     RESTART = 'restart.mat';               % Init data from restart data

    VIZ = 'yes';                            % Dynamically vizualize results

%     STOP_WITH = 'S';            % stopping criterion is simplices
    STOP_WITH = 'V';            % stopping criterion is vertices (samples)
    
    % Way to refine simplexes
    REFINE_TYPE = 'UNIFORM';    % based on volume of simplex only
%     REFINE_TYPE = 'PWC_LES';    % error estimate, large edge split
%     REFINE_TYPE = 'PWC_MP';     % error estimate, midpoint split
%     REFINE_TYPE = 'PWC_LES_RAND';    % error estimate, large edge split
                                
    ERROR_LIMIT = 1e-12;         % default value for error estimate if is zero
    p = 1;                      % local polynomial interpolant degree
%-------------------------------------------------------------------------%

    if exist('RESTART') == 0
        
        % Declare important parameters:
        N = 2;                      % Dimensionality of stochastic space
        % Test Function (distribution)
%         TEST_FUNC = 'NORMAL';       % normalized guassian in all dimensions  
%         TEST_FUNC = 'COS_2';        % cos^2 in all dimensions
%         TEST_FUNC = 'ND_HUMP';      % two normal distributions in each dim
%         TEST_FUNC = 'SINC';         % matlab sinc function
%         TEST_FUNC = 'STEP_ISH';     % 0.1 for x<pi/10 in all dimensions
        TEST_FUNC = 'J_ARCTAN';       % Jeroen arctan problem
%         TEST_FUNC = 'BETA';       % Jeroen arctan problem

        % Generate initial distribution of points:
        % corners of every dimension (2^N points)
        NODE = zeros(2^N,N);
        NODE = de2bi( (0:(2^N-1))' );

        % Add midpoint sample (for now)
        NODE(2^N + 1,:) = mean(NODE);

        % Build simplex array
        % each row is the i_th simplex, values listed are the indices
        % of the NODEs which define the simplex
        SIMPLEX = delaunayn(NODE);
        
        known_points = 2^N + 1;
        
        switch STOP_WITH
            case 'S'
                disp(cat(2,'Initialized # of simplexes: ',num2str(size(SIMPLEX,1))));
                n_stop = input('Stop after how many total simplexes? ');
                if isempty(n_stop)
                    n_stop=100;
                end
            case 'V'
                disp(cat(2,'Initialized # of vertices: ',num2str(size(NODE,1))));
                n_stop = input('Stop after how many total vertices (samples)? ');
                if isempty(n_stop)
                    n_stop=100;
                end
            otherwise
                disp('Unrecognized stopping condition');
        end
    
        TOT_ERR = zeros(n_stop - 2*factorial(N),1);

        % Evaluate solution and error estimate
        SOLUTION = zeros(size(NODE,1),1);

        % loop over all vertices
        for v=1:size(NODE,1)
            % evaulate solution
            SOLUTION(v) = exact_solution(NODE(v,:));
        end



        % Build initial error estimate for simplices
        % piecewise constant, value is equal to error in refinement which
        % generated that simplex
        ERROR = zeros(size(SIMPLEX,1),1);
        ERR_HAT = zeros(size(SIMPLEX,1),1);
        for s=1:size(SIMPLEX,1)
            ERROR(s) = 10e6;
            
            ERR_HAT(s) = ERROR(s) ^ ((p+1)/N);
        end
        
        % while loop counter initialization
        counter=1;
        
        % plot init
        old_stop = 1;
        
    elseif exist('RESTART') == 1
        % start from data specified in RESTART
        load(RESTART);
        clear DATA;
        
        % for plotting purposes - old error
        old_stop=size(TOT_ERR,1);
        
        N = size(NODE,2);
        
        switch STOP_WITH
            case 'S'
                disp(cat(2,'Initialized # of simplexes: ',num2str(size(SIMPLEX,1))));
                n_stop = input('Stop after how many total simplexes? ');
                if isempty(n_stop)
                    n_stop=size(SIMPLEX,1);
                end
            case 'V'
                disp(cat(2,'Initialized # of vertices: ',num2str(size(NODE,1))));
                n_stop = input('Stop after how many total vertices (samples)? ');
                if isempty(n_stop)
                    n_stop=size(NODE,1);
                end
            otherwise
                disp('Unrecognized stopping condition');
        end
    
        TOT_ERR = cat(1,TOT_ERR,zeros(n_stop-2*factorial(N)-old_stop,1));
    end
    
    % Compute volume of all simplices in stochastic space
    VOL = zeros(size(SIMPLEX,1),1);
    % loop over all simplices
    for s=1:size(SIMPLEX,1) 
        VOL(s) = volume_of_simplex(s,SIMPLEX,N,NODE);
    end
    
    STOP = 0;

    while STOP == 0
        
        % original number of nodes/samples
        n_node_old = size(NODE,1);
        
        switch REFINE_TYPE
            case 'UNIFORM'
                
                % retrieve index of simplex with largest volume
                [value,index] = max(VOL);
             
                % refine that simplex
                [SIMPLEX,NODE,ERROR,ERR_HAT] = ...
                    refine_simplex_LES(index,SIMPLEX,N,NODE,ERROR,ERR_HAT);
                
            case 'PWC_MP'
                % New point added is in the middle of the current simplex
                % Refinement based on piecewise constant error estimate
                SIMPLEX_ERROR = zeros(size(SIMPLEX,1),1);
                % loop over all simplices
                for s=1:size(SIMPLEX,1) 
                    SIMPLEX_ERROR(s) = max(ERR_HAT(s) * VOL(s),ERROR_LIMIT);
                end
                
                % retrieve index of simplex with largest volume
                [value,index] = max(SIMPLEX_ERROR);       
                
                % refine that simplex
                [SIMPLEX,NODE,ERROR,ERR_HAT] = ...
                    refine_simplex_MP(index,SIMPLEX,N,NODE,ERROR,ERR_HAT);
                
            case 'PWC_LES'
                % Refinement based on piecewise constant error estimate
                SIMPLEX_ERROR = zeros(size(SIMPLEX,1),1);
                % loop over all simplices
                for s=1:size(SIMPLEX,1) 
                    SIMPLEX_ERROR(s) = max(ERR_HAT(s) * VOL(s),ERROR_LIMIT);
                end
                
                % retrieve index of simplex with largest volume
                [value,index] = max(SIMPLEX_ERROR);   
                
                % refine that simplex
                [SIMPLEX,NODE,ERROR,ERR_HAT] = ...
                    refine_simplex_LES(index,SIMPLEX,N,NODE,ERROR,ERR_HAT);
                
            case 'PWC_LES_RAND'
                % Refinement based on piecewise constant error estimate
                SIMPLEX_ERROR = zeros(size(SIMPLEX,1),1);
                % loop over all simplices
                for s=1:size(SIMPLEX,1) 
                    SIMPLEX_ERROR(s) = max(ERR_HAT(s) * VOL(s),ERROR_LIMIT);
                end
                
                total_error = sum(SIMPLEX_ERROR);
                
                % normalize error distribution, randomly sample discrete
                % space
                index = discretesample(SIMPLEX_ERROR/total_error,1);      
                
                % refine that simplex
                [SIMPLEX,NODE,ERROR,ERR_HAT] = ...
                    refine_simplex_LES(index,SIMPLEX,N,NODE,ERROR,ERR_HAT);
                
            otherwise
                disp('REFINE_TYPE not supported');
        end
        
        % Re-compute volume of all simplices in stochastic space
        VOL = zeros(size(SIMPLEX,1),1);
        % loop over all simplices
        for s=1:size(SIMPLEX,1) 
            VOL(s) = volume_of_simplex(s,SIMPLEX,N,NODE);
        end
        TOT_VOL = sum(VOL);

        % update total error estimate for this refinement loop
        if counter > 2*factorial(N)
            TOT_ERR(counter-2*factorial(N)) = sum(ERR_HAT .* (VOL ./ TOT_VOL));
        end
        
        % If new point was added, generate SOLUTION and update total error
        if n_node_old == (size(NODE,1)-1)
            SOLUTION(size(NODE,1)) = exact_solution(NODE(size(NODE,1),:));

            SAMPLE_ERR(size(NODE,1)-(2^N + 1),1) = sum(ERR_HAT .* (VOL ./ TOT_VOL));
            disp(cat(2,'Sample #  : ',num2str(size(NODE,1))));
        end
        
        disp(cat(2,'Refinement #  : ',num2str(counter)));
            
%         STOP = input('Keep refining? (1 to STOP):');
%         if isempty(STOP)
%             STOP=0;
%         end
        
        if N <= 3 && N > 0 && exist('VIZ') == 1
            % Plot Graphics Page
            plot_figs(NODE,SIMPLEX,N,SOLUTION,TOT_ERR,old_stop);
        end
        
        % stopping criterion
        switch STOP_WITH
            case 'S'
                if size(SIMPLEX,1) >= n_stop
                    STOP = 1;
                end
            case 'V'
                if size(NODE,1) >= n_stop
                    STOP = 1;
                end
        end

        counter = counter + 1;
    end

    if N <= 3 && N > 0
        % Plot Graphics Page
        plot_figs(NODE,SIMPLEX,N,SOLUTION,TOT_ERR,old_stop);
    end
    
    disp(cat(2,'# of Refinements: ',num2str(counter-1)));
    disp(cat(2,'# of Simplices  : ',num2str(size(SIMPLEX,1))));
    disp(cat(2,'# of Vertices   : ',num2str(size(NODE,1))));
    disp(cat(2,'Total Vol       : ',num2str(TOT_VOL)));
    disp(cat(2,'Total Error     : ',num2str(sum(ERR_HAT .* (VOL./TOT_VOL)))));
    
    %------------------------------------
    % Error Estimate After Refinement
    %------------------------------------

    % Compute RMS error
    rms = sqrt( sum( (VOL./ TOT_VOL) .* ERR_HAT.^2 ) );
    disp(cat(2,'RMS Error       : ',num2str(rms)));
    
    % Build data structure to hold results
    DATA.NODE=NODE;
    DATA.SIMPLEX=SIMPLEX;
    DATA.ERROR=ERROR;
    DATA.ERR_HAT=ERR_HAT;
    DATA.SOLUTION=SOLUTION;
    DATA.TEST_FUNC=TEST_FUNC;
    DATA.TOT_ERR=TOT_ERR;
    DATA.counter=counter;
    DATA.SAMPLE_ERR=SAMPLE_ERR;
    
    if N <= 3 && N > 0
        % Plot Graphics Page
        plot_convergence(N,SAMPLE_ERR);
    end
    
    save restart.mat -struct DATA
end

%% Functions


%--------------------------------
% Interpolate solution
%--------------------------------
function [value] = interp_solution(x,i,SIMPLEX,N,NODE)
    % interpolate solution at location x based on nodes which make up
    % simplex i
    X = [];
    Y = [];
    for n=1:size(SIMPLEX,2)
       X = cat(1,X,NODE(SIMPLEX(i,n),:)); 
       Y = cat(1,Y,exact_solution(NODE(SIMPLEX(i,n),:)));
    end
            
    switch N
        case 1
            value = interp1(X,Y,x);
            
        case 2
            
            F = TriScatteredInterp(X,Y);
            value = F(x); 
        case 3

            F = TriScatteredInterp(X,Y);
            value = F(x);
            
        otherwise
            disp('Interp in ndim > 3 not yet supported');
    end

end

%--------------------------------
% Evaluate Exact Solution at
% given sample point
%--------------------------------
function [solution] = exact_solution(x)
    % pass in x which is location in n-D
    global TEST_FUNC;

    switch TEST_FUNC
        case 'NORMAL'
            mean = 0.5;
            std = 0.15;

            % maximum of normal distribution, to nomralize result to [0,1]
            max = normpdf(mean,mean,std);

            % for each dimension of the node's location, compute normal
            % distribution and generate normalized product
            solution = 1;
            for j=1:size(x,2)
                solution = solution * normpdf(x(j),mean,std) / max;
            end

        case 'ND_HUMP'
            mean1 = 0.25;
            mean2 = 0.75;
            std1 = 0.1;
            std2 = 0.2;

            % maximum of normal distribution, to normalize result to [0,1]
            max = normpdf(mean1,mean1,std1)+normpdf(mean1,mean2,std2);

            % for each dimension of the node's location, compute normal
            % distribution and generate normalized product
            solution = 1;
            for j=1:size(x,2)
                solution = solution * ( normpdf(x(j),mean1,std1) + normpdf(x(j),mean2,std2) )...
                    / max;
            end
        case 'COS_2'
            % for each dimension of the node's location, compute normal
            % distribution and generate normalized product
            solution = 1;
            for j=1:size(x,2)
                solution = solution * ( cos(x(j)*pi)^2 );
            end
            
        case 'SINC'
            solution = 1;
            for j=1:size(x,2)
                solution = solution * sinc(2*pi*(x(j)));
            end
            
        case 'STEP_ISH'
            if sum(x < (0.1*pi)) == size(x,2)
                % in all dimensions is less than 0.5
                solution = 0.1;
            else
                solution = 1.0;
            end
            
        case 'J_ARCTAN'
            x_star = ones(1,size(x,2));
            
            solution = atan(dot(x,x_star) + x_star(1)^2);
            
        case 'BETA'
            solution = 1;
            for j=1:size(x,2)
                solution = solution * betapdf(x(j),4,4);
            end
            
        otherwise
            disp('TEST_FUNC not supported');
    end
end

%--------------------------------
% Refine simplex using 
% Largest Edge Splitting
%--------------------------------
function [new_SIMPLEX,NODE,new_ERROR,new_ERR_HAT] = ...
    refine_simplex_MP(i,SIMPLEX,N,NODE,ERROR,ERR_HAT)    
    % given a refinement simplex (i), generate updated SIMPLEX and NODE
    % lists by adding a point at the midpoint of the simplex
    global TEST_FUNC;
    global p;
    
    % new point is the average of all nodes in the current simplex
    new_point = 0;
    for n=1:size(SIMPLEX,2)
        new_point = new_point + NODE(SIMPLEX(i,n),:);
    end
    new_point = new_point / size(SIMPLEX,2);
    
    % search to make sure point does not already exist in NODE array
    [n_value, n_index] = max(ismember(NODE,new_point,'rows'));
    
    if n_value == 1
        % Yes, does already exist
        new_node = n_index;
    else
        % No, does not already exist
        % add point to NODE array
        new_node = size(NODE,1)+1;
        NODE(new_node,:) = new_point(1,:);
    end
    
    % reference simplex volume
    old_VOL = volume_of_simplex(i,SIMPLEX,N,NODE);
    
    % perform delaunay triangulation on JUST the refined element
    % collect nodes/vertices for old simplex
    simplex_points = [];
    for j=1:size(SIMPLEX,2)
        simplex_points = cat(1,simplex_points,NODE(SIMPLEX(i,j),:));
    end
    % add new vertex
    simplex_points = cat(1,simplex_points,NODE(new_node,:));
    SIMPLEX_temp = delaunayn(simplex_points,{'QJ','QbB'});
    
    % temp indexing of new points
    index = cat(2,SIMPLEX(i,:),new_node);
    % generate SIMPLEX list of new simplices using NODE index, instead of 
    % local cell index
    SIMPLEX_new = zeros(size(SIMPLEX_temp,1),size(SIMPLEX_temp,2));
    ERROR_new = zeros(size(SIMPLEX_temp,1),1);
    ERR_HAT_new = zeros(size(SIMPLEX_temp,1),1);
    for r=1:size(SIMPLEX_temp,1)
       for c=1:size(SIMPLEX_temp,2)
           SIMPLEX_new(r,c) = index(SIMPLEX_temp(r,c));
       end

       % get volume of new simplex
       new_VOL = volume_of_simplex(r,SIMPLEX_temp,N,NODE);
       
       % Compute error between estimate and exact
       ERROR_new(r) = abs(interp_solution(new_point(1,:),i,SIMPLEX,N,NODE)...
        - exact_solution(NODE(new_node,:)));
       
       ERR_HAT_new(r) = ERROR_new(r) * ( (new_VOL/old_VOL)^((p+1)/N) );
    end      

    new_SIMPLEX = [];
    new_ERROR = [];
    new_ERR_HAT = [];
    % Pop current row, append new tesselation
    for row=1:size(SIMPLEX,1)
       if row ~= i
          new_SIMPLEX = cat(1,new_SIMPLEX,SIMPLEX(row,:));
          new_ERROR = cat(1,new_ERROR,ERROR(row));
          new_ERR_HAT = cat(1,new_ERR_HAT,ERR_HAT(row));
       end

    end

    new_ERROR = cat(1,new_ERROR,ERROR_new);
    new_ERR_HAT = cat(1,new_ERR_HAT,ERR_HAT_new);
    new_SIMPLEX = cat(1,new_SIMPLEX,SIMPLEX_new);
end


%--------------------------------
% Refine simplex using 
% Largest Edge Splitting
%--------------------------------
function [new_SIMPLEX,NODE,new_ERROR,new_ERR_HAT] = ...
    refine_simplex_LES(i,SIMPLEX,N,NODE,ERROR,ERR_HAT)    
    % given a refinement simplex (i), generate updated SIMPLEX and NODE
    % lists
    global TEST_FUNC;
    global p;
    
    % N choose 2 edges for a given simplex
    edges = combntns(1:(N+1),2);
    
    % Compute all edge lengths
    edge_length = zeros(size(edges,1),1);
    
    for e = 1:size(edges,1)
       edge_length(e) = sqrt( sum( (NODE(SIMPLEX(i,edges(e,2)),:)-NODE(SIMPLEX(i,edges(e,1)),:)).^2 ));
    end
    
    % retrieve index of longest edge
    [le_value,le_index] = max(edge_length);
    
    new_point = 0.5 * (NODE(SIMPLEX(i,edges(le_index,2)),:) + NODE(SIMPLEX(i,edges(le_index,1)),:));
    
    % search to make sure point does not already exist in NODE array
    [n_value, n_index] = max(ismember(NODE,new_point,'rows'));
    
    if n_value == 1
        % Yes, does already exist
        new_node = n_index;
    else
        % No, does not already exist
        % add point to NODE array
        new_node = size(NODE,1)+1;
        NODE(new_node,:) = new_point(1,:);
    end
    
    % reference simplex volume
    old_VOL = volume_of_simplex(i,SIMPLEX,N,NODE);
    
    % perform delaunay triangulation on JUST the refined element
    % collect nodes/vertices for old simplex
    simplex_points = [];
    for j=1:size(SIMPLEX,2)
        simplex_points = cat(1,simplex_points,NODE(SIMPLEX(i,j),:));
    end
    % add new vertex
    simplex_points = cat(1,simplex_points,NODE(new_node,:));
    SIMPLEX_temp = delaunayn(simplex_points,{'QJ','QbB'});
    
    % temp indexing of new points
    index = cat(2,SIMPLEX(i,:),new_node);
    % generate SIMPLEX list of new simplices using NODE index, instead of 
    % local cell index
    SIMPLEX_new = zeros(size(SIMPLEX_temp,1),size(SIMPLEX_temp,2));
    ERROR_new = zeros(size(SIMPLEX_temp,1),1);
    ERR_HAT_new = zeros(size(SIMPLEX_temp,1),1);
    for r=1:size(SIMPLEX_temp,1)
       for c=1:size(SIMPLEX_temp,2)
           SIMPLEX_new(r,c) = index(SIMPLEX_temp(r,c));
       end

       % get volume of new simplex
       new_VOL = volume_of_simplex(r,SIMPLEX_new,N,NODE);
       
       % Compute error between estimate and exact
       ERROR_new(r) = abs(interp_solution(new_point(1,:),i,SIMPLEX,N,NODE)...
        - exact_solution(NODE(new_node,:)));
       
       ERR_HAT_new(r) = ERROR_new(r) * ( (new_VOL/old_VOL)^((p+1)/N) );
    end      

    new_SIMPLEX = [];
    new_ERROR = [];
    new_ERR_HAT = [];
    % Pop current row, append new tesselation
    for row=1:size(SIMPLEX,1)
       if row ~= i
          new_SIMPLEX = cat(1,new_SIMPLEX,SIMPLEX(row,:));
          new_ERROR = cat(1,new_ERROR,ERROR(row));
          new_ERR_HAT = cat(1,new_ERR_HAT,ERR_HAT(row));
       end

    end

    new_ERROR = cat(1,new_ERROR,ERROR_new);
    new_ERR_HAT = cat(1,new_ERR_HAT,ERR_HAT_new);
    new_SIMPLEX = cat(1,new_SIMPLEX,SIMPLEX_new);
end

%--------------------------------
% Compute Volume of simplex
%--------------------------------
function [volume] = volume_of_simplex(i,SIMPLEX,N,NODE)
    % Corresponds to Eq. (6) in SSC Paper
    vertices=SIMPLEX(i,:);
    
    temp=[];
    for j=2:size(SIMPLEX,2)
        temp = cat(1,temp,NODE(vertices(j-1),:)-NODE(vertices(j),:));
    end

    volume = abs(det(temp))/factorial(N);
end

%--------------------------------
% Plot Figures/Graphs
%--------------------------------
function plot_figs(NODE,SIMPLEX,N,SOLUTION,TOT_ERR,old_stop)
    
    % Init plot figure
    figure(1);
    
    % Scatter plot
    switch N
        case 1        
            % SSC refinement
            subplot(1,2,1);hold on;cla;
            [spoints,sI]=sort(NODE);
            surf(...
                cat(2,spoints(:,1),spoints(:,1)),...                                % X
                cat(2,-0.1*ones(length(spoints),1),0.1*ones(length(spoints),1)),... % Y
                cat(2,SOLUTION(sI),SOLUTION(sI)));                                  % Z
            axis equal;view(20,30);
            hold off;
            
            % Total Error Plot
            subplot(1,2,2);hold on;cla;
            plot(1:old_stop,TOT_ERR(1:old_stop),'-r','LineWidth',4);
            plot(old_stop:(size(TOT_ERR,1)-2*factorial(N)),...
                TOT_ERR(old_stop:(size(TOT_ERR,1)-2*factorial(N))),'-b','LineWidth',4);
            hold off;
            
        case 2
            % Scatter plot
            subplot(1,3,1);hold on;cla;
            scatter(NODE(:,1),NODE(:,2),'ok','filled')
            axis equal;axis([-0.1 1.1 -0.1 1.1]);
            hold off;
            
            % SSC refinement
            subplot(1,3,2);hold on;cla;
            trisurf(SIMPLEX,NODE(:,1),NODE(:,2),SOLUTION);
            axis equal;axis([-0.1 1.1 -0.1 1.1]);view(20,30);
            hold off;
            
            % Total Error Plot
            subplot(1,3,3);hold on;cla;
            plot(1:old_stop,TOT_ERR(1:old_stop),'-r','LineWidth',4);
            plot(old_stop:(size(TOT_ERR,1)-2*factorial(N)),...
                TOT_ERR(old_stop:(size(TOT_ERR,1)-2*factorial(N))),'-b','LineWidth',4);
            hold off;
            
        case 3
            % Scatter plot
            subplot(1,3,1);hold on;cla;
            scatter3(NODE(:,1),NODE(:,2),NODE(:,3),'ok','filled')
            axis equal;axis([-0.1 1.1 -0.1 1.1 -0.1 1.1]);view(20,30);
            grid on;
            hold off;
            
            % SSC refinement
            subplot(1,3,2);hold on;cla;
            tetramesh(SIMPLEX,NODE,'FaceAlpha',0.3);
            axis equal;axis([-0.1 1.1 -0.1 1.1 -0.1 1.1]);view(20,30);
            hold off;
            
            % Total Error Plot
            subplot(1,3,3);hold on;cla;
            plot(1:old_stop,TOT_ERR(1:old_stop),'-r','LineWidth',4);
            plot(old_stop:(size(TOT_ERR,1)-2*factorial(N)),...
                TOT_ERR(old_stop:(size(TOT_ERR,1)-2*factorial(N))),'-b','LineWidth',4);
            hold off;
        otherwise
            disp('Dimension not supported');     
    end
    
%     figure(2);
%     hold on;cla;
%     loglog(1:old_stop,TOT_ERR(1:old_stop),'-r','LineWidth',4);
%     loglog(old_stop:(size(TOT_ERR,1)-2*factorial(N)),...
%         TOT_ERR(old_stop:(size(TOT_ERR,1)-2*factorial(N))),'-b','LineWidth',4);
%     hold off;
end

%--------------------------------
% Plot Figures/Graphs
%--------------------------------
function plot_convergence(N,SAMPLE_ERR)
    
    % Init plot figure
 
    stop = size(SAMPLE_ERR,1);
    
    figure(2);
    loglog((2^N+1):stop,SAMPLE_ERR((2^N+1):stop),'-k','LineWidth',4);
    hold on;
    % plot slope = -1. -2 lines
    loglog([1 stop]+2^N,SAMPLE_ERR(2^N+1)*[1 10^(-2*log10(stop))],'--r','LineWidth',2); 
    loglog([1 stop]+2^N,SAMPLE_ERR(2^N+1)*[1 10^(-1*log10(stop))],'--b','LineWidth',2); 
    grid on;
    hold off;
    
end

function x = discretesample(p, n)
% Samples from a discrete distribution
%
%   x = discretesample(p, n)
%       independently draws n samples (with replacement) from the 
%       distribution specified by p, where p is a probability array 
%       whose elements sum to 1.
%
%       Suppose the sample space comprises K distinct objects, then
%       p should be an array with K elements. In the output, x(i) = k
%       means that the k-th object is drawn at the i-th trial.
%       
%   Remarks
%   -------
%       - This function is mainly for efficient sampling in non-uniform 
%         distribution, which can be either parametric or non-parametric.         
%
%       - The function is implemented based on histc, which has been 
%         highly optimized by mathworks. The basic idea is to divide
%         the range [0, 1] into K bins, with the length of each bin 
%         proportional to the probability mass. And then, n values are
%         drawn from a uniform distribution in [0, 1], and the bins that
%         these values fall into are picked as results.
%
%       - This function can also be employed for continuous distribution
%         in 1D/2D dimensional space, where the distribution can be
%         effectively discretized.
%
%       - This function can also be useful for sampling from distributions
%         which can be considered as weighted sum of "modes". 
%         In this type of applications, you can first randomly choose 
%         a mode, and then sample from that mode. The process of choosing
%         a mode according to the weights can be accomplished with this
%         function.
%
%   Examples
%   --------
%       % sample from a uniform distribution for K objects.
%       p = ones(1, K) / K;
%       x = discretesample(p, n);
%
%       % sample from a non-uniform distribution given by user
%       x = discretesample([0.6 0.3 0.1], n);
%
%       % sample from a parametric discrete distribution with
%       % probability mass function given by f.
%       p = f(1:K);
%       x = discretesample(p, n);
%

%   Created by Dahua Lin, On Oct 27, 2008
%

%% parse and verify input arguments

assert(isfloat(p), 'discretesample:invalidarg', ...
    'p should be an array with floating-point value type.');

assert(isnumeric(n) && isscalar(n) && n >= 0 && n == fix(n), ...
    'discretesample:invalidarg', ...
    'n should be a nonnegative integer scalar.');

%% main

% process p if necessary

K = numel(p);
if ~isequal(size(p), [1, K])
    p = reshape(p, [1, K]);
end

% construct the bins

edges = [0, cumsum(p)];
s = edges(end);
if abs(s - 1) > eps
    edges = edges * (1 / s);
end

% draw bins

rv = rand(1, n);
c = histc(rv, edges);
ce = c(end);
c = c(1:end-1);
c(end) = c(end) + ce;

% extract samples

xv = find(c);

if numel(xv) == n  % each value is sampled at most once
    x = xv;
else                % some values are sampled more than once
    xc = c(xv);
    d = zeros(1, n);
    dv = [xv(1), diff(xv)];
    dp = [1, 1 + cumsum(xc(1:end-1))];
    d(dp) = dv;
    x = cumsum(d);
end

% randomly permute the sample's order
x = x(randperm(n));
end
