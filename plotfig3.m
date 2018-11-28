fig=figure();
hold on
a={'576.mat','5760.0001.mat','5760.00001.mat'};
label={'0.001','0.0001','0.00001'};

d=load('5760.00001.mat');
for i=1:3
    d=load(a{i})
    cdfplot(d.save_data(:,5))
end
legend(label)
% title('5 minutes, MEM call number CDF')
xlabel('Number of ReMEM call')
saveas(fig,'figs/5min_memcall_cdf_all.png')
% saveas(fig,'figs/5min_memcall.pdf')

