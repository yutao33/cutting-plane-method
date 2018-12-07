function [best_a,i] = volumetric_center_method(C, d, varargin)
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



% Some debug and info variable
feasible_found = 0;
num_till_feasible = 0; %number of iterations until first feasible solution
num_inside = 0; %number of times point inside ellipsoid is found after first feasible solution
num_outside = 0; %number of times point outside ellipsoid is used after first feasible solution

initial_R=1;
P_A = -[-eye(n);eye(n)];
P_b = -[zeros(n,1);ones(n,1)*initial_R];
a=ones(n,1)*0.5*initial_R;

error_distance=initial_R;

best_a = a;
best_obj = 0;

%% 
delta=1e-4;
varepsilon=1e-3*delta;


fprintf('N=%d\n',N);
for i=1:N

    [m,~]=size(P_A); % TODO check n;
    z=a;
    sigmma=zeros(1,m);
    for ii=1:m
        sigmma(ii)=sigmma_i(ii,P_A,P_b,z);
    end
    [min_sigmma, min_index]=min(sigmma);
    if min_sigmma>=varepsilon  % case 1: add plane
        
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
        % [c, gamma, result] = separation_oracle_omniscient(C, d, a);
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
    
    tmp=-c';
    tmplength=norm(tmp);
    tmp=tmp./tmplength;
    
    gamma = tmp*z-max(1e-7,sqrt((tmp/H(P_A,P_b,z)*tmp')/2*sqrt(delta*varepsilon)));
    [lia, loc]=ismember(tmp,P_A,'rows');
    if lia
        if P_b(loc)<gamma
            P_b(loc)=gamma;
        end
    else
        P_A=[P_A;tmp];
        P_b=[P_b;gamma];
    end
    z=Newton_method(P_A,P_b,z,varepsilon);
    if any(isnan(z))
        error('isnan')
    end
    % gradient_F(P_A,P_b,z)
    error_distance=max(abs(a-z));
    if error_distance<eps
       disp(error_distance) 
    end
    a=z;
    
    % error_distance=2*m*(det(H(P_A,P_b,z)))^(-1/2/n)
    
    else  % case 2: remote plane min_sigmma
    
    P_A(min_index,:)=[];
    P_b(min_index,:)=[];
    z=Newton_method(P_A,P_b,z,varepsilon);
    if any(isnan(z))
        error('isnan')
    end
    % gradient_F(P_A,P_b,z)
    % error_distance=sum(abs(a-z))
%     if error_distance<eps
%        disp(error_distance) 
%     end
    a=z;
    % error_distance=2*m*(det(H(P_A,P_b,z)))^(-1/2/n)
    
    end

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


function r=Newton_method(A,b,z,varepsilon)
    for j=1:ceil(30*log(2*varepsilon^(-4.5)))
%         e=0.18*(Q(A,b,z)\gradient_F(A,b,z));
        [Q,gF]=Q_gF(A,b,z);
        e=0.18*(Q\gF);
        if any(isnan(e))
            error('e isnan')
        end
        z=z-e;
        if all(abs(e)<=1e-10)
            j;
            break
        end
    end
    r=z;
end

function [Q,gF]=Q_gF(A,b,x)
    [m,n]=size(A);
    Q=zeros(n,n);
    gF=zeros(n,1);
    % Hv=H(A,b,x);
    Hv=zeros(n,n);
    tmpv=zeros(m,1);
    tmpv_2=zeros(m,1);
    for i=1:m
        ai=A(i,:)';
        bi=b(i);
        tmp1 = ai'*x-bi;
        tmp1_2 = tmp1^2;
        tmpv(i)=tmp1;
        tmpv_2(i)=tmp1_2;
        Hv=Hv+ai*ai'./tmp1_2;
    end
    for i=1:m
        ai=A(i,:)';
        % bi=b(i);
        tmp1 = tmpv(i);
        tmp1_2 = tmpv_2(i);
        sigmma_i=(ai'/Hv*ai)./tmp1_2;
        tmp = sigmma_i * ai;
        Q=Q + tmp*ai'./tmp1_2;
        gF=gF - tmp./tmp1;
    end
end

function r=gradient_F(A,b,x)
    [m,n]=size(A);
    r=zeros(n,1);
    for i=1:m
        ai=A(i,:)';
        r=r-sigmma_i(i,A,b,x)*ai./(ai'*x-b(i));
    end
end

function r=Q(A,b,x)
    [m,n]=size(A);
    r=zeros(n,n);
    for i=1:m
        ai=A(i,:)';
        bi=b(i);
        r=r+sigmma_i(i,A,b,x)*(ai*ai')./(ai'*x-bi)^2;
    end
end

function r=sigmma_i(i,A,b,x)
    ai=A(i,:)';
    bi=b(i);
    r=(ai'/H(A,b,x)*ai)./(ai'*x-bi)^2;
end

function r=H(A,b,x)
    [m,n]=size(A);
    r=zeros(n,n);
    for i=1:m
        ai=A(i,:)';
        bi=b(i);
        
        r=r+ai*ai'./(ai'*x-bi)^2;
    end
end