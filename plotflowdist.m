data_dir='D:\MyWork\ellipsoid-data\constraints-intervals';

for num=[1] % 288 192 144 96 48]
    name=num2str(num);
    subdir=strcat(data_dir,'/',name,'/*.mat');
    allfile = struct2cell(dir(subdir));
    [k len]=size(allfile);
    save_data=[];
    filename_list={};
    for i=1:len
        filename=allfile{1,i};
        filename = strcat(data_dir,'/',name,'/', filename);
        [result, cow,row]=test_case(filename);
        if result
            save_data=[save_data; [cow,row]];
            filename_list{end+1}=filename;
        end
    end
    save('result/flow_dist.mat','save_data','filename_list');
    v=save_data(:,1)';
    v1=min(v);
    v2=200 %max(v);
    x=v1:v2;
    y=zeros(1,v2-v1+1);
    for vv=v
        if ismember(vv,x)
            y(vv-v1+1)=y(vv-v1+1)+1;
        end
    end
    fig=figure();
    bar(x,y);
end

function [result, cow,row]=test_case(filename)
disp(filename);
data=csvread(filename);
cow = size(data,2)-1;
row = size(data,1);
result = true;
end