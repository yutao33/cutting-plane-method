global special_count_1
special_count_1=[];
global MEM_count MEM_count1 MEM_count2
MEM_count=0;
MEM_count1=0;
MEM_count2=0;

global eps_global
eps_global=0.001;
o = [1 0.9];
C = [1 0;
    0 1;
    -1 0;
    0 -1;
    1 1];
d = [1;1;0;0;1.5];
t1=clock();
[a,iter] = randomwalk_method(C,d,o,'deep','optimize',1,'radius',2,'plot_fig',0,'eps',eps_global);
time=etime(clock(),t1);

count1=MEM_count1;
count2=MEM_count2;