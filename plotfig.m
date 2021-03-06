% 2880-config8 eps=0.01 r=0.05 R=2
% 2880-config9 eps=0.0001 r=0.001 R=2

if ~exist('savefig','var')
   savefig=false; 
end

% set(0,'DefaultFigureVisible', 'off')
% name='2880-config8';
name='2880-18-11-30-09-24-45-config';
filename=strcat('result/',name,'.mat');

fig=figure();
d=load(filename);
cdfplot(d.save_data(:,7))
title('')
xlabel('Number of ReSEP calls')
ylabel('')
set(gcf,'position',[403   397   366   269])
axis normal
if savefig
    saveas(fig,strcat('figs/',name,'-sep-cdf.pdf'))
end


fig=figure();
row=size(d.save_data,1);
nn=13;
a=ones(row,nn)*NaN;
for i=1:row
    n=d.save_data(i,1);
    a(i,n)=d.save_data(i,7);
end
boxplot(a,[1:nn])
xlabel('Number of flows')
ylabel('Number of ReSEP calls')
set(gcf,'position',[403   397   366   269])
if savefig
    saveas(fig,strcat('figs/',name,'-sep-box.pdf'))
end



fig=figure();
d=load(filename);
cdfplot(d.save_data(:,5))
title('')
xlabel('Number of total ReMEM calls')
ylabel('')
set(gcf,'position',[403   397   366   269])
if savefig
    saveas(fig,strcat('figs/',name,'-mem-cdf.pdf'))
end

fig=figure();
row=size(d.save_data,1);
% nn=7;
a=ones(row,nn)*NaN;
for i=1:row
    n=d.save_data(i,1);
    a(i,n)=d.save_data(i,5);
end
boxplot(a,[1:nn])
xlabel('Number of flows')
ylabel('Number of total ReMEM calls')
set(gcf,'position',[403   397   366   269])
if savefig
    saveas(fig,strcat('figs/',name,'-mem-box.pdf'))
end



latency = d.save_data(:,7) .* (35+rand(row,1)*70);

fig=figure();
cdfplot(latency)
title('')
xlabel('Latency (ms)')
ylabel('')
set(gcf,'position',[403   397   366   269])
if savefig
    saveas(fig,strcat('figs/',name,'-latency-cdf.pdf'))
end


fig=figure();
row=size(d.save_data,1);
% nn=7;
a=ones(row,nn)*NaN;
for i=1:row
    n=d.save_data(i,1);
    a(i,n)=latency(i);
end
boxplot(a,[1:nn])
xlabel('Number of flows')
ylabel('Latency (ms)')
set(gcf,'position',[403   397   366   269])
if savefig
    saveas(fig,strcat('figs/',name,'-latency-box.pdf'))
end


constrainsts=d.save_data(:,2);

% d=load(filename);
% cdfplot(constrainsts)
% title('')
% xlabel('Number of Constraints')
% ylabel('')
% set(gcf,'position',[403   397   366   269])
% saveas(fig,strcat('figs/',name,'-constraints-cdf.pdf'))

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


