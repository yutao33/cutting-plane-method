name='ss';
fig=figure();

% d=load(name);
cdfplot(special_count_1(:,2))
title('')
xlabel(' x label ')
ylabel('')
set(gcf,'position',[403   397   366   269])
axis normal
saveas(fig,strcat('figs/',name,'-cdf.pdf'))
% saveas(fig,'figs/5min_memcall.pdf')

fig=figure();
row=size(special_count_1,1);
nn=7;
a=ones(row,nn)*NaN;
for i=1:row
    n=special_count_1(i,1);
    a(i,n)=special_count_1(i,2);
end
boxplot(a,[1:nn])
xlabel('xlable')
ylabel('ylabel')
set(gcf,'position',[403   397   366   269])
saveas(fig,strcat('figs/',name,'-box.pdf'))