config=struct();
config.savefig=false;
config.savefigpath='figs2/';
config.figsize=[403   397   400   300];

% filename_list={'1-18-12-02-16-39-21-config-ellipsoid','1-18-12-02-15-23-23-config-rand','1-18-12-07-17-17-05-config-vo-1-15'};
filename_list={'1-18-12-07-22-00-10-config-el-1-15','1-18-12-07-20-45-23-config-ra-1-15','1-18-12-07-17-17-05-config-vo-1-15'};
savename_list={'ellipsoid','random-walk','volumetric-center'};

for i=1:3
    filename=filename_list{i};
    config.name=savename_list{i};
    fullfilename=strcat('result/',filename,'.mat');
    data = load(fullfilename);

    cow_range=1:15;
    sepnum=data.save_data(:,7);
    cow=data.save_data(:,1);
    select_matrix=ismember(cow,cow_range);
    select_sepnum = sepnum(select_matrix);
    
    %CDF_sepnum_one(select_sepnum, config);
    %BOX_sepnum_one(cow, sepnum, cow_range, config);
    
    memnum=data.save_data(:,5);
    select_memnum = memnum(select_matrix);
    %CDF_memnum_one(select_memnum, config);
    BOX_memnum_one(cow, memnum, cow_range, config);
    %latency_plot(cow, sepnum, cow_range, config);
end
% close all;

%%
function CDF_sepnum_one(data, config)
fig=figure();
cdfplot(data)
title('')
xlabel('Number of ReSEP calls')
ylabel('')
set(gcf,'position',config.figsize)
axis normal
if config.savefig
    saveas(fig,strcat(config.savefigpath,config.name,'-sepnum-cdf.pdf'))
end
end

%%
function BOX_sepnum_one(cow, sepnum, cow_range, config)
fig=figure();
nn=length(cow_range);
a=[];
for i=1:length(cow)
    cow_num=cow(i);
    [loc,loi]=ismember(cow_num,cow_range);
    if loc
        aa=ones(1,nn)*NaN;
        aa(loi)=sepnum(i);
        a=[a;aa];
    end
end
boxplot(a,cow_range);
xlabel('Number of flows');
ylabel('Number of ReSEP calls');
set(gcf,'position',config.figsize);
if config.savefig
    saveas(fig,strcat(config.savefigpath,config.name,'-sep-box.pdf'))
end
end

%%
function CDF_memnum_one(data, config)
fig=figure();
cdfplot(data);
title('')
xlabel('Number of total ReMEM calls')
ylabel('')
set(gcf,'position',config.figsize);
if config.savefig
    saveas(fig,strcat(config.savefigpath,config.name,'-mem-cdf.pdf'))
end
end

%%
function BOX_memnum_one(cow, memnum, cow_range, config)
fig=figure();
nn=length(cow_range);
a=[];
for i=1:length(cow)
    cow_num=cow(i);
    [loc,loi]=ismember(cow_num,cow_range);
    if loc
        aa=ones(1,nn)*NaN;
        aa(loi)=memnum(i);
        a=[a;aa];
    end
end
boxplot(a,cow_range);
xlabel('Number of flows')
ylabel('Number of total ReMEM calls')
set(gcf,'position',config.figsize)
if config.savefig
    saveas(fig,strcat(config.savefigpath,config.name,'-mem-box.pdf'))
end
end

%%
function latency_plot(cow, sepnum, cow_range, config)
row=length(cow);
latency = sepnum .* (35+rand(row,1)*70);

select_matrix=ismember(cow,cow_range);
latency_selected=latency(select_matrix);

fig=figure();
cdfplot(latency_selected)
title('')
xlabel('Latency (ms)')
ylabel('')
set(gcf,'position',config.figsize)
if config.savefig
    saveas(fig,strcat(config.savefigpath,config.name,'-latency-cdf.pdf'))
end


fig=figure();
nn=length(cow_range);
a=[];
for i=1:row
    cow_num=cow(i);
    [loc,loi]=ismember(cow_num,cow_range);
    if loc
        aa=ones(1,nn)*NaN;
        aa(loi)=latency(i);
        a=[a;aa];
    end
end

boxplot(a,cow_range)
xlabel('Number of flows')
ylabel('Latency (ms)')
set(gcf,'position',config.figsize)
if config.savefig
    saveas(fig,strcat(config.savefigpath,config.name,'-latency-box.pdf'))
end
end


%%
function constrainsts_num()

constrainsts=d.save_data(:,2);


fig=figure();
minn=0 %min(constrainsts);
maxn=max(constrainsts);
delta=10;
group=ceil((maxn+0.00001)/delta);
a=ones(row,group)*NaN;
for i=1:row
    n=ceil((constrainsts(i,1)-minn)/delta+0.00001);
    if n==0
        n=1;
    elseif n==group+1
        n=group;
    end
    a(i,n)=d.save_data(i,7);
end
a
% num = [minn+1/2*delta:delta:maxn]
num=cell(group,1);
for i=1:group
    num(i)={sprintf('%d~%d',minn+(i-1)*delta,minn+i*delta)};
end
num
boxplot(a,num)
xtb=get(gca,'XTickLabel');
xt=get(gca,'XTick')
yt=get(gca,'YTick')
ytextp=yt(1)*ones(1,length(xt)); 
text(xt+0.1,ytextp-50*0.5,xtb,'HorizontalAlignment','right','rotation',45,'fontsize',9); % 50
set(gca,'xticklabel','')
xla=xlabel('Number of Constraints');
yla=ylabel('Number of ReSEP calls');
set(gcf,'position',[403   397   366   269])
xla_po=get(xla,'Position')-[0 85*0.5 0]; % 85
set(xla,'Position',xla_po);
if savefig
    saveas(fig,strcat('figs/',name,'-constraints-sepcall-box.pdf'))
end

%%
fig=figure();
a=ones(row,group)*NaN;
for i=1:row
    n=ceil((constrainsts(i,1)-minn)/delta+0.00001);
    if n==0
        n=1;
    elseif n==group+1
        n=group;
    end
    a(i,n)=d.save_data(i,5);
end
boxplot(a, num)
xtb=get(gca,'XTickLabel');
xt=get(gca,'XTick')
yt=get(gca,'YTick')
ytextp=yt(1)*ones(1,length(xt)); 
text(xt+0.1,ytextp-10000*0.6,xtb,'HorizontalAlignment','right','rotation',45,'fontsize',9);
set(gca,'xticklabel','')
xla=xlabel('Number of Constraints');
yla=ylabel('Number of total ReMEM calls');
set(gcf,'position',[403   397   366   269])
xla_po=get(xla,'Position')-[0 8500*1.2 0]
set(xla,'Position',xla_po);
if savefig
    saveas(fig,strcat('figs/',name,'-constraints-memcall-box.pdf'))
end

end


