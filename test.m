% test script
clear
clc
nowdatestr = datestr(now,'yy-mm-dd-HH-MM-SS');
diary(strcat('log/log-',nowdatestr,'.txt'))
diary on;

data_dir='D:\MyWork\ellipsoid-data\constraints-intervals';

global looplimit
looplimit=zeros(1,30);

for num=[1] % 288 192 144 96 48]
    name=num2str(num);
    subdir=strcat(data_dir,'\',name,'\*.mat');
    allfile = struct2cell(dir(subdir));
    [k len]=size(allfile);
    save_data=[];
    count1_s=[];
    count2_s=[];
    
    global special_count_1
    special_count_1=[];
    
    for i=1:len
        filename=allfile{1,i};
        filename = strcat(data_dir,'\',name,'\', filename);
        [result, cow,row,linprog_best,ellipsoid_best, MEM_call, time, count1, count2]=test_case(filename);
        if result
            save_data=[save_data; [cow,row,linprog_best,ellipsoid_best, MEM_call, time]];
            count1_s=[count1_s, count1];
            count2_s=[count2_s, count2];
        end
    end
    save_data(:,7)=count1_s';
    save_data(:,8)=count2_s';
    label={'cow','row','linprog_best','ellipsoid_best', 'MEM_call', 'time','count1_s','count2_s'};
    config={'eps=0.01 r=0.01 R=2 1-7'};
    save(strcat('result/',name,'-',nowdatestr,'-config.mat'),'save_data','label','config');
end
diary off




function [result, cow,row,linprog_best,ellipsoid_best, MEM_call, time, count1, count2]=test_case(filename)
result=false;
cow=-1;
row=-1;
linprog_best=-1;
ellipsoid_best=-1;
MEM_call=-1;
time=-1;
count1=-1;
count2=-1;

disp(filename);
data=csvread(filename);

cow = size(data,2);

if cow>15 || cow<=9
    return
end

% global looplimit
% if looplimit(cow)>4
%     return
% end
% looplimit(cow)=looplimit(cow)+1;

C=data(:,1:cow-1);
d=data(:,cow);
d=d/10000000;
if d>1
    error('error range');
end

row=size(C,1);
cow=cow-1;
C(row+1:row+cow,:)=-eye(cow);
d(row+1:row+cow,:)=zeros(cow,1);
o=zeros(1,cow)+1;

options = optimoptions('linprog','Display', 'off','Algorithm','dual-simplex','OptimalityTolerance',1e-8,'ConstraintTolerance',1e-8,'MaxIterations',1e6);
[x,fval,exitflag,output]=linprog(-o,C,d,[],[],[],[],[],options);
linprog_best=-fval;
fprintf('linprog best=%f\n',linprog_best);
fprintf('exitflag=%d\n',exitflag)
if exitflag~=1 || linprog_best>100000
    C
    d
    o
    disp(output)
    disp('linprog failed');
    return
end

global MEM_count MEM_count1 MEM_count2
MEM_count=0;
MEM_count1=0;
MEM_count2=0;

global eps_global
eps_global=0.0001;
% deep cut
t1=clock();
[a,iter] = volumetric_center_method(C,d,o,'deep','optimize',1,'radius',2,'plot_fig',0,'eps',eps_global);
time=etime(clock(),t1);

count1=MEM_count1;
count2=MEM_count2;

format longG
fprintf('Best feasible point after iteration %i:\n', iter);
disp(a)
ellipsoid_best=o*a;
fprintf('Objective Value: %f\n', ellipsoid_best);
fprintf('MEM_count=%d\n',MEM_count);
fprintf('MEM_count1=%d\n',MEM_count1);
fprintf('MEM_count2=%d\n',MEM_count2);
fprintf('cow=%d\n',cow);
MEM_call=count1+count2;
% error('d')
result=true;
end