function [result, c,gamma] = separation_oracle_2(C, d, a)
% if all(C*a <= d)
%     c = a; % a is in polytope, so return a
%     gamma = 0; % gamma not needed
% else
%     % find any violated constraint
%     j = find(C*a > d, 1);
%     % return violated inequality as separating hyperplane
%     c = C(j,:)';
%     % value of inequality (needed for deep cut)
%     gamma = d(j);
% end
global C_g d_g
C_g=C;
d_g=d;
if isnan(a)
    error('a isnan')
end
[result,g,const]=weaksep(a);

c=g;

gamma=c'*a+const;

if ~all(isreal(c))
   error('c is not real') 
end

% if result
%     c=g;
%     gamma=c'*a;
% else
%     c=g;
%     gamma=c'*a+const;
% end
end

function test=MEM(a)
global C_g d_g MEM_count
global special_count
test=all(C_g*a<=d_g);
MEM_count=MEM_count+1;
special_count=special_count+1;
end

function [result,g,const]=weaksep(x)
n=length(x);

global eps_global
epsilon=eps_global; %0.0001;
rou=epsilon; % TODO check it

r=0.001;  % TODO check it
R=2; % sqrt(n); 

kappa=R/r;

% rou=sqrt(n^(7/6)*kappa^2*)
rou=1;

global MEM_count1
MEM_count1=MEM_count1+1;

global special_count
special_count=0;

if MEM(x)
    result=true;
    g=x;
    const=0;
elseif x'*x>R
    result=false;
    g=x;
    const=0;
else
    r1=n^(1/6) * epsilon^(1/3) * R^(2/3) / kappa;
    p=@(d) h_func(x,d);
    g=approxsubgradient(p,zeros(n,1),r1,4*epsilon, 3*kappa);
    while all(g==0)
        disp('repeat')
        g=approxsubgradient(p,zeros(n,1),r1,4*epsilon, 3*kappa);
    end
    result=false;
    const1=50/rou * n^(7/6) * R^(2/3) * kappa * epsilon^(1/3);
%     fprintf('const = %f\n',const1)
    const=0;
    global special_count_1    
    special_count_1=[special_count_1; [n,special_count]];
end





end

function h=h_func(x, d)
global MEM_count2

high=1.1;   % NOTE: special for x>=0
low=0;
if ~MEM(d)
    error('d not in MEM')
end
if MEM(x)
    error('x in MEM')
end
while high-low>0.000001
    MEM_count2=MEM_count2+1;
    
    mid=(high+low)/2;
    if MEM(d+mid*x)
        low=mid;
    else
        high=mid;
    end
end
mid=(high+low)/2;
h=-mid*sqrt(x'*x);
end

function g_tidle = approxsubgradient(p, x, r1, epsilon, L)
% L=1000; % TODO checkit
n=length(x);
r2 = sqrt(epsilon*r1/(sqrt(n)*L));
% random_select=@(x,r)();
y=random_select(x,r1);
z=random_select(y,r2);
g_tidle=zeros(n,1);
for i=1:n
    si=z;
    si(i)=y(i)-r2;
    di=z;
    di(i)=y(i)+r2;
    g_tidle(i)=(p(di)-p(si))/(2*r2);
end

if all(g_tidle==0)
%     error('dd')
    disp('dd')
end
end

function y = random_select(x,r)
y=x+r*(rand(length(x),1)*0.8+0.2);
end
