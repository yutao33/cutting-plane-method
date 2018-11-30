function [best_a,i] = randomwalk_method(C, d, varargin)
%% Parse input parameters
p = inputParser;

defaultCut = 'central';
validCuts = {'central','deep','shallow'};
checkCut = @(x) any(validatestring(x,validCuts));

addRequired(p,'C',@ismatrix);
addRequired(p,'d',@isvector);
addOptional(p,'objective',@isvector);
addOptional(p,'cut',defaultCut,checkCut);

addParameter(p,'optimize',0,@isnumeric);
addParameter(p,'eps',1e-5,@isnumeric);
addParameter(p,'numiter',0,@isnumeric);
addParameter(p,'radius',0,@isnumeric);
addParameter(p,'ignore_blowup',0,@isnumeric);
addParameter(p,'center',[0;0],@isvector);

addParameter(p,'plot_fig',0,@isnumeric);
addParameter(p,'plot_iter',[0 1],@isvector);
addParameter(p,'plot_title','',@ischar);
addParameter(p,'plot_separating',1,@isnumeric);
addParameter(p,'plot_gradient',1,@isnumeric);

addParameter(p,'plot_pause',0,@isnumeric);

parse(p,C,d,varargin{:});

% Assign input parameters
o = p.Results.objective;
eps = p.Results.eps;
OPTIMIZE = p.Results.optimize;
ignore_blowup = p.Results.ignore_blowup;

plot_fig = p.Results.plot_fig;
plot_iter = p.Results.plot_iter;
plot_title = p.Results.plot_title;
plot_sep = p.Results.plot_separating;
plot_grd = p.Results.plot_gradient;
plot_pause = p.Results.plot_pause;

FONTSIZE = 11;
[~,n] = size(C);

% Use input parameters if set or calculate default ones
if p.Results.numiter == 0
    %N = 2*n*((2*n+1)*input_size(C)+n*input_size(d)-n^3); %according to Grötschel 3.1.32
    N = 50 * (n+1)^2 * input_size(C,d); %according to Grötschel 3.1.37
    %N = ceil(10 * n^2 * (2 * log2(n) + 5 * input_size(C,d))); %according to Korte p. 99
    %N = 5 * n^2 * ceil(log( (6*R^2*max( norm(o), 1 )) / (r*eps) )); %according to Korte p. 101
else
    N = p.Results.numiter;
end
if p.Results.radius == 0
    R = max( n * 2^(2*(input_size(C,d)-n^2)) , 1); %according to Grötschel 3.1.32 & p. 80 avoid sqrt!
    %R = max( n * (1 + 2^(4 * input_size(C,d))) , 1); %according to Korte p. 99
else
    R = p.Results.radius^2;
end
%precision = 8*N; % OPEN: precision parameter?

% Initialization of ellipsoid matrix and center
A = R * eye(n); % take R^2 above to avoid calc of sqrt in line 32 (see Grötschel p. 80)
if size(p.Results.center,1) ~= n
    a = zeros(n,1);
else
    a = p.Results.center;
end

% Save best found feasible solution
best_a = a;
best_obj = 0;

% Some debug and info variable
feasible_found = 0;
num_till_feasible = 0; %number of iterations until first feasible solution
num_inside = 0; %number of times point inside ellipsoid is found after first feasible solution
num_outside = 0; %number of times point outside ellipsoid is used after first feasible solution

initial_R=1;
randomwalk_area_A = [-eye(n);eye(n)];
randomwalk_area_b = [zeros(n,1);ones(n,1)*initial_R];

randomwalk_result=cprnd(n*8+1,randomwalk_area_A,randomwalk_area_b);
a=mean(randomwalk_result)';
% a=ones(n,1)*0.5*initial_R;
error_distance=initial_R;

%% 

fprintf('N=%d\n',N)
for i=1:N
    if ~all(a>=0)
        n=length(a);
        result=false;
        j=find(a<0,1);
        c=zeros(n,1);
        c(j)=-1;
        gamma=0;
    elseif feasible_found==1 && o*a<best_obj
        % patch
        result=false;
        c=-o';
        gamma=-best_obj;
    else
        [result, c, gamma] = separation_oracle_2(C, d, a);
%        [c, gamma] = separation_oracle_omniscient(C, d, a);
%         if all(a==c)
%             result=true;
%         else
%             result=false;
%         end
    end
    % Oracle just returned input variable "a" (all constraints are satisfied)
%     if all(a == c)
    if result
        % First feasible solution found
        if feasible_found == 0
            num_till_feasible = i-1; %solution was already found in previous iteration
            feasible_found = 1; %set to 1 at first found feasible solution
        end
		% Optimization via sliding objective method
        if OPTIMIZE
            % Add objective as constraint for sliding objective optimization
            c = -o'; % o' for min // -o' for max
            if error_distance > eps
                gamma = -o*a; % o for min // -o for max
                num_inside = num_inside + 1; %for info purpose: count iterations inside polytope
                if o*a > best_obj
                    best_a = a;
                    best_obj = o*a;
                end
            else
                best_a = a;
                break;
            end
        else
            best_a = a;
            break;
        end
    else
        if feasible_found == 1
           num_outside = num_outside + 1; 
        end
        if OPTIMIZE
            if error_distance > eps

            else
                best_a = a;
                break;
            end
        else

        end
    end
    
    % add c*x< gamma
    tmp=c';
    tmplength=norm(tmp);
    tmp=tmp./tmplength;
    gamma=gamma/tmplength;
    [lia, loc]=ismember(tmp,randomwalk_area_A,'rows');
    if lia
        if randomwalk_area_b(loc)>gamma
            randomwalk_area_b(loc)=gamma;
        end
    else
        randomwalk_area_A=[randomwalk_area_A;tmp];
        randomwalk_area_b=[randomwalk_area_b;gamma];        
    end
    half_result=randomwalk_result(randomwalk_result*tmp'<=gamma,:);
    half_center=mean(half_result,1);
    options=struct();
    options.x0=half_center';
    randomwalk_result=cprnd(n*8+1,randomwalk_area_A,randomwalk_area_b,options);
    for ii=1:(n*8+1)
        if any(randomwalk_area_A*randomwalk_result(ii,:)'>randomwalk_area_b+1e-8)
            error('randomwalk failed');
        end
    end
    
    a=mean(randomwalk_result)';
    error_distance=norm(max(randomwalk_result)-min(randomwalk_result));

end

%% Result of the ellipsoid method
% Print some information about the output and plot the final result.

% Check if solution in last iteration is feasible (in case the very last iteration found a feasible solution)
if feasible_found == 0
    [c,~] = separation_oracle_omniscient(C, d, a);
    if all(a == c)
        feasible_found = 1;
        num_till_feasible = i;
        best_a = a;
    end
end

% Output info
if feasible_found == 1
    fprintf('First feasible solution found in iteration %i!\n', num_till_feasible);
    if OPTIMIZE
        fprintf('Optimal value after iteration %i with %i iterations outside polytope and %i iterations inside polytope.\n', i-1, num_outside, num_inside);
    end
else
    fprintf('Warning: Solution not feasible after final iteration %i!\n', i);
end

end