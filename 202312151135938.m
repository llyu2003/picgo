%{
llyu@04.24   
新思路：计算不同降水等级下各月能见度样本中位数形成的91站*12个月
以此为基础，分别计算年和季节，得到91站*5年四季，季节(年)值=对应月份想加除以3(12)
区域阈值的方法：nanmean(91站*12个月)  和   nanmean(91站*5年四季)
%}
%% 1.cal-hourly降水分成小中大暴下各站medianVIS--各月
%VIS自动观测始于2016
%%tag-yyyy-mm-dd-hh-t(摄氏度)-pre（mm）-rh(%)-vis(m)
%仔细观察数据，时刻从21时到20时，除了降水为逐时数据，其他都是一天四次
%观测的四个时刻为02/08/14/20，显然为BJT
%思路：先等级样本，再取中位数
clear,clc
pathin='Z:\data\obs\2001-2021pre-t-rh\pro\时间调整后\';
pathout='Z:\data\obs\2001-2021pre-t-rh\pro\';
matpath='Z:\data\obs\2001-2021pre-t-rh\pro\matfile\';
pltpath='Z:\data\obs\2001-2021pre-t-rh\fig\';
load([matpath,'stn91_classify.mat']);
%包含变量station_info，station_types， 'plain_stns','hill_stns','mountain_stns'
xstation_types={'平原站','丘陵站','高山站'};%对应{'plain','hill','mountain'}

PREclass.str{1,1}='逐时降水划分为小雨、中雨、大雨和暴雨四个等级';
PREclass.str{2,1}='降水等级标准-国标GB/T25892-2012';
PREclass.cn={'无雨','小雨','中雨','大雨','暴雨'};
PREclass.eng={'noRain','light','median','heavy','storm'};
PREclass.below=[0,0.1,1.6,7.0,15.0];
PREclass.up=[0.09,1.5,6.9,14.9,999.9];

ntypes=size(station_types,2);
nclass=size(PREclass.below,2);
for ic=1:nclass
for itype=1:ntypes
    tagsinfo=eval([station_types{1,itype},'_stns']); %tag-lon-lat-ht
    tags=tagsinfo(:,1);
    xthres_vis=zeros(numel(tags),16)*nan;%tag-lon-lat-ht-minvis[1-12month]
    for id=1:numel(tags)
        raw=dlmread([pathin,num2str(tags(id)),'hourly-T-PRE-RH-VIS.txt']);
        %%tag-yyyy-mm-dd-hh-t(摄氏度)【6】-pre（mm）【7】-rh(%)【8】-vis(m)【9】
        %VIS自动观测始于2016,截取1月1日0时起  气象时调整后
        line_2016=find(raw(:,2)==2016 & raw(:,3)==1 & raw(:,4)==1 & raw(:,5)==0); 
        line_end=find(raw(:,2)==2020 & raw(:,3)==12 & raw(:,4)==31 & raw(:,5)==23);
        sample2x=raw(line_2016:line_end,:);%vis没有取2021年，因为2021年数据用于检验了        
        %sample2=sample2x(sample2x(:,7)>=0.1,:);%研究vis，待分级降水的时序样本,后来发现可以不用它

        vec_vis=zeros(1,12)*nan;%minvis[1to12]          
        for im=1:12  %月份
            %先取等级样本---依据降水量
            class_sample=sample2x((sample2x(:,7)>=PREclass.below(ic) & sample2x(:,7)<=PREclass.up(ic)),:);
            %tag-yyyy-mm-dd-hh-t(摄氏度)【6】-pre（mm）【7】-rh(%)【8】-vis(m)【9】            
            
            %{
            %再去异常值---三倍标准差
            xvec=class_sample(:,9); %待处理的样本
            xave=nanmean(xvec);   xstd=nanstd(xvec); %均值和均方根
            up_val=xave+3*xstd;  below_val=xave-3*xstd; %上限和下限
            inx=find(xvec<=up_val & xvec>=below_val); %正常值对应的索引
            tmpx=class_sample(inx,:); % 正常值
            iny=find(xvec>up_val | xvec<below_val); %异常值对应的索引
            tmpyy=class_sample(iny,:); % to check 异常值  
            %最后取逐月
            sample_for_vis=tmpx(tmpx(:,3)==im,:);%研究vis，所有降水的时序样本            
            %}            
            %最后取逐月
            sample_for_vis=class_sample(class_sample(:,3)==im,:);%研究vis，所有降水的时序样本
            %==============================================
            %-----------求某等级降水某月对应的vis最小值-----------            
            if ~isempty(sample_for_vis)
                %minvis=nanmin(sample_for_vis(:,9));%出现降水对应的vis阈值
                minvis=nanmedian(sample_for_vis(:,9));%出现降水对应的vis阈值
            else
                minvis=nan;
            end
            vec_vis(1,im)=round(minvis,0);   %vec_vis(1,im)=minvis;

        end %逐月求阈值 im
        xthres_vis(id,:)=[tagsinfo(id,:),vec_vis];
        %xthres_vis一共16列：tag-lon-lat-ht-minvis[1-12月]
    end %某类台站循环  id
    eval(['minvis_',num2str(itype),'=xthres_vis;']);%tag-lon-lat-ht-minvis[1-12月]
end %台站分类 itype
thres_vis=[minvis_1;minvis_2;minvis_3];%三类台站的结果垂直叠加，没采用追加形式
str1={'tag','lon','lat','ht','1月','2月','3月','4月','5月','6月','7月','8月','9月','10月','11月','12月'};

%str1={'tag','lon','lat','ht','年','春','夏','秋','冬'};
%str1={'tag','lon','lat','ht','春','夏','秋','冬'};
%dlmwrite([matpath,'各站四季vis阈值.txt'],str1{:},'delimiter','\t');
%'precision','%.4f',

fname2=strcat('中位数from',num2str(PREclass.below(ic)),'to',num2str(PREclass.up(ic)),'mm各站逐月vis阈值2016-2020.txt');
fid2=fopen([matpath,fname2],'w+');%w+是重写内容
for ii=1:size(str1,2)
fprintf(fid2,'%s\t',str1{1,ii});%'%s\t'表明输出字符串后加一个tab键，好跟后面的字符串分开
end
fprintf(fid2,'\n');
fclose(fid2);
dlmwrite([matpath,fname2],thres_vis,'delimiter','\t','-append');
%-append表示末尾追加，在第一行str1换行后追加thres_vis（thres_vis已经是所有台站结果了）
%～～～～～～～～～～～～～～～～～～～～～～··
end %降水等级 ic

%% 2.cal-hourly降水分成6个等级各站medianVIS--各月
%VIS自动观测始于2016
%%tag-yyyy-mm-dd-hh-t(摄氏度)-pre（mm）-rh(%)-vis(m)
%仔细观察数据，时刻从21时到20时，除了降水为逐时数据，其他都是一天四次
%观测的四个时刻为02/08/14/20，显然为BJT
%出现降水时刻对应的VIS历史最低值
%思路：先等级样本，再取中位数

clear,clc
pathin='Z:\data\obs\2001-2021pre-t-rh\pro\时间调整后\';
pathout='Z:\data\obs\2001-2021pre-t-rh\pro\';
matpath='Z:\data\obs\2001-2021pre-t-rh\pro\matfile\';
pltpath='Z:\data\obs\2001-2021pre-t-rh\fig\';
load([matpath,'stn91_classify.mat']);
%包含变量station_info，station_types， 'plain_stns','hill_stns','mountain_stns'
xstation_types={'平原站','丘陵站','高山站'};%对应{'plain','hill','mountain'}

PREclass.str{1,1}='逐时降水划分为六个等级';
PREclass.str{2,1}='降水等级标准-自定义';
PREclass.cn={'[0.1,0.5)','[0.5,2)','[2,5)','[5,10)','[10,20)','[20,999.9)'};
PREclass.eng={'lev1','lev2','lev3','lev4','lev5','lev6'};
PREclass.below=[0.1,0.5,2,5,10,20];
PREclass.up=[0.4,1.9,4.9,9.9,19.9,999.9];

ntypes=size(station_types,2);
nclass=size(PREclass.below,2);
for ic=1:nclass
for itype=1:ntypes
    tagsinfo=eval([station_types{1,itype},'_stns']); %tag-lon-lat-ht
    tags=tagsinfo(:,1);
    xthres_vis=zeros(numel(tags),16)*nan;%tag-lon-lat-ht-minvis[1-12month]
    for id=1:numel(tags)
        raw=dlmread([pathin,num2str(tags(id)),'hourly-T-PRE-RH-VIS.txt']);
        %%tag-yyyy-mm-dd-hh-t(摄氏度)【6】-pre（mm）【7】-rh(%)【8】-vis(m)【9】
        %VIS自动观测始于2016,截取1月1日0时起  气象时调整后
        line_2016=find(raw(:,2)==2016 & raw(:,3)==1 & raw(:,4)==1 & raw(:,5)==0); 
        line_end=find(raw(:,2)==2020 & raw(:,3)==12 & raw(:,4)==31 & raw(:,5)==23);
        sample2x=raw(line_2016:line_end,:);%vis没有取2021年，因为2021年数据用于检验了        
        %sample2=sample2x(sample2x(:,7)>=0.1,:);%研究vis，待分级降水的时序样本,后来发现可以不用它

        vec_vis=zeros(1,12)*nan;%minvis[1to12]          
        for im=1:12  %月份
            %先取等级样本---依据降水量
            class_sample=sample2x((sample2x(:,7)>=PREclass.below(ic) & sample2x(:,7)<=PREclass.up(ic)),:);
            %tag-yyyy-mm-dd-hh-t(摄氏度)【6】-pre（mm）【7】-rh(%)【8】-vis(m)【9】            
            
            %{
            %再去异常值---三倍标准差
            xvec=class_sample(:,9); %待处理的样本
            xave=nanmean(xvec);   xstd=nanstd(xvec); %均值和均方根
            up_val=xave+3*xstd;  below_val=xave-3*xstd; %上限和下限
            inx=find(xvec<=up_val & xvec>=below_val); %正常值对应的索引
            tmpx=class_sample(inx,:); % 正常值
            iny=find(xvec>up_val | xvec<below_val); %异常值对应的索引
            tmpyy=class_sample(iny,:); % to check 异常值             
            %最后取逐月
            sample_for_vis=tmpx(tmpx(:,3)==im,:);%研究vis，所有降水的时序样本
            %}
            
            sample_for_vis=class_sample(class_sample(:,3)==im,:);%研究vis，所有降水的时序样本
            %==============================================
            %-----------求某等级降水某月对应的vis最小值-----------            
            if ~isempty(sample_for_vis)
                minvis=nanmedian(sample_for_vis(:,9));%出现降水对应的vis阈值
                %minvis=nanmin(sample_for_vis(:,9));%出现降水对应的vis阈值
            else
                minvis=nan;
            end
            vec_vis(1,im)=round(minvis,0);   %vec_vis(1,im)=minvis;

        end %逐月求阈值 im
        xthres_vis(id,:)=[tagsinfo(id,:),vec_vis];
        %xthres_vis一共16列：tag-lon-lat-ht-minvis[1-12月]
    end %某类台站循环  id
    eval(['minvis_',num2str(itype),'=xthres_vis;']);%tag-lon-lat-ht-minvis[1-12月]
end %台站分类 itype
thres_vis=[minvis_1;minvis_2;minvis_3];%三类台站的结果垂直叠加，没采用追加形式
str1={'tag','lon','lat','ht','1月','2月','3月','4月','5月','6月','7月','8月','9月','10月','11月','12月'};

%str1={'tag','lon','lat','ht','年','春','夏','秋','冬'};
%str1={'tag','lon','lat','ht','春','夏','秋','冬'};
%dlmwrite([matpath,'各站四季vis阈值.txt'],str1{:},'delimiter','\t');
%'precision','%.4f',

fname2=strcat('中位数from',num2str(PREclass.below(ic)),'to',num2str(PREclass.up(ic)),'mm各站逐月vis阈值2016-2020.txt');
fid2=fopen([matpath,fname2],'w+');%w+是重写内容
for ii=1:size(str1,2)
fprintf(fid2,'%s\t',str1{1,ii});%'%s\t'表明输出字符串后加一个tab键，好跟后面的字符串分开
end
fprintf(fid2,'\n');
fclose(fid2);
dlmwrite([matpath,fname2],thres_vis,'delimiter','\t','-append');
%-append表示末尾追加，在第一行str1换行后追加thres_vis（thres_vis已经是所有台站结果了）
%～～～～～～～～～～～～～～～～～～～～～～··
end %降水等级 ic

%% 3.cal-hourly降水不分等级各站medianVIS--各月
%VIS自动观测始于2016
%%tag-yyyy-mm-dd-hh-t(摄氏度)-pre（mm）-rh(%)-vis(m)
%仔细观察数据，时刻从21时到20时，除了降水为逐时数据，其他都是一天四次
%观测的四个时刻为02/08/14/20，显然为BJT
%出现降水时刻对应的VIS历史最低值
%思路：先等级样本，再取中位数
clear,clc
pathin='Z:\data\obs\2001-2021pre-t-rh\pro\时间调整后\';
pathout='Z:\data\obs\2001-2021pre-t-rh\pro\';
matpath='Z:\data\obs\2001-2021pre-t-rh\pro\matfile\';
pltpath='Z:\data\obs\2001-2021pre-t-rh\fig\';
load([matpath,'stn91_classify.mat']);
%包含变量station_info，station_types， 'plain_stns','hill_stns','mountain_stns'
xstation_types={'平原站','丘陵站','高山站'};%对应{'plain','hill','mountain'}

PREclass.str{1,1}='不分等级';
PREclass.str{2,1}='大于等于0.1或0.5mm';
PREclass.cn={'大于等于0.1','大于等于0.5'};
PREclass.eng={'ge0.1','ge0.5'};
PREclass.below=[0.1,0.5];
PREclass.up=[999.9,999.9];

ntypes=size(station_types,2);
nclass=size(PREclass.below,2);
for ic=1:nclass
for itype=1:ntypes
    tagsinfo=eval([station_types{1,itype},'_stns']); %tag-lon-lat-ht
    tags=tagsinfo(:,1);
    xthres_vis=zeros(numel(tags),16)*nan;%tag-lon-lat-ht-minvis[1-12month]
    for id=1:numel(tags)
        raw=dlmread([pathin,num2str(tags(id)),'hourly-T-PRE-RH-VIS.txt']);
        %%tag-yyyy-mm-dd-hh-t(摄氏度)【6】-pre（mm）【7】-rh(%)【8】-vis(m)【9】
        %VIS自动观测始于2016,截取1月1日0时起  气象时调整后
        line_2016=find(raw(:,2)==2016 & raw(:,3)==1 & raw(:,4)==1 & raw(:,5)==0); 
        line_end=find(raw(:,2)==2020 & raw(:,3)==12 & raw(:,4)==31 & raw(:,5)==23);
        sample2x=raw(line_2016:line_end,:);%vis没有取2021年，因为2021年数据用于检验了        
        %sample2=sample2x(sample2x(:,7)>=0.1,:);%研究vis，待分级降水的时序样本,后来发现可以不用它

        vec_vis=zeros(1,12)*nan;%minvis[1to12]          
        for im=1:12  %月份
            %先取等级样本---依据降水量
            class_sample=sample2x((sample2x(:,7)>=PREclass.below(ic) & sample2x(:,7)<=PREclass.up(ic)),:);
            %tag-yyyy-mm-dd-hh-t(摄氏度)【6】-pre（mm）【7】-rh(%)【8】-vis(m)【9】            
            
            %{
            %再去异常值---三倍标准差
            xvec=class_sample(:,9); %待处理的样本
            xave=nanmean(xvec);   xstd=nanstd(xvec); %均值和均方根
            up_val=xave+3*xstd;  below_val=xave-3*xstd; %上限和下限
            inx=find(xvec<=up_val & xvec>=below_val); %正常值对应的索引
            tmpx=class_sample(inx,:); % 正常值
            iny=find(xvec>up_val | xvec<below_val); %异常值对应的索引
            tmpyy=class_sample(iny,:); % to check 异常值              
            %最后取逐月
            sample_for_vis=tmpx(tmpx(:,3)==im,:);%研究vis，所有降水的时序样本
            %}
            
            sample_for_vis=class_sample(class_sample(:,3)==im,:);%研究vis，所有降水的时序样本
            %==============================================
            %-----------求某等级降水某月对应的vis最小值-----------                
            if ~isempty(sample_for_vis)
                minvis=nanmedian(sample_for_vis(:,9));%出现降水对应的vis阈值
                %minvis=nanmin(sample_for_vis(:,9));%出现降水对应的vis阈值
            else
                minvis=nan;
            end
            vec_vis(1,im)=round(minvis,0);   %vec_vis(1,im)=minvis;

        end %逐月求阈值 im
        xthres_vis(id,:)=[tagsinfo(id,:),vec_vis];
        %xthres_vis一共16列：tag-lon-lat-ht-minvis[1-12月]
    end %某类台站循环  id
    eval(['minvis_',num2str(itype),'=xthres_vis;']);%tag-lon-lat-ht-minvis[1-12月]
end %台站分类 itype
thres_vis=[minvis_1;minvis_2;minvis_3]; %三类台站的结果垂直叠加，没采用追加形式
str1={'tag','lon','lat','ht','1月','2月','3月','4月','5月','6月','7月','8月','9月','10月','11月','12月'};

%str1={'tag','lon','lat','ht','年','春','夏','秋','冬'};
%str1={'tag','lon','lat','ht','春','夏','秋','冬'};
%dlmwrite([matpath,'各站四季vis阈值.txt'],str1{:},'delimiter','\t');
%'precision','%.4f',

fname2=strcat('中位数from',num2str(PREclass.below(ic)),'to',num2str(PREclass.up(ic)),'mm各站逐月vis阈值2016-2020.txt');
fid2=fopen([matpath,fname2],'w+');%w+是重写内容
for ii=1:size(str1,2)
fprintf(fid2,'%s\t',str1{1,ii});%'%s\t'表明输出字符串后加一个tab键，好跟后面的字符串分开
end
fprintf(fid2,'\n');
fclose(fid2);
dlmwrite([matpath,fname2],thres_vis,'delimiter','\t','-append');
%-append表示末尾追加，在第一行str1换行后追加thres_vis（thres_vis已经是所有台站结果了）
%～～～～～～～～～～～～～～～～～～～～～～··
end %降水等级 ic  

%% 4.计算不同条件下各站年和四季的阈值，基于91台站*12个月medianVIS--年四季 
%基于区域数据：91台站*12个月。计算冬季，某站12月，1月，2月，算数平均作为该站冬季的阈值，
%其他季节亦如此，如此就有91台站*四季；再以此为依据，某季节的91个阈值做算术平均作为该季节区域阈值
clear,clc
list=dir('Z:\data\obs\2001-2021pre-t-rh\pro\matfile\2023*');
%手动点击list看看，有name对应三个文件夹，folder为上级路径{'Z:\data\obs\2001-2021pre-t-rh\pro\matfile'}
nf=size(list,1); %获取文件夹的个数  
matpath=list(1).folder;  

for i=1:nf
    %subfoldername=list.name{i,:};%wrong wrong wrong
    subfoldername=list(i).name;
    subpath=strcat(matpath,'\',subfoldername,'\');
     sublist2=dir([subpath,'中位数*逐月*.txt']); 
    sort_nat_name2=sort_nat({sublist2.name}); %使用sort_nat函数对文件排序
    %=======
    %   seasonal
    %=======
    for ix=1:size(sublist2,1)  %txt文件的个数 
        a=dlmread([subpath,sort_nat_name2{ix}],'\t',1,0);%用排序后的文件名读入
        %a=dlmread([subpath,sublist1(ix).name],'\t',1,0);%左上角的起始坐标是R=0,C=0 ，得到91*9的double数据  
        %a2=importdata([subpath,sublist1(ix).name],'\t',1);%1表示取第一行以下，得到一个结构数据
        %a is tag-lat-lon-ht-1到12个月
        
        %===cal ==============
        ann=nanmean(a(:,5:end),2);
        spring=nanmean(a(:,7:9),2); 
        summer=nanmean(a(:,10:12),2); 
        autumn=nanmean(a(:,13:15),2); 
        winter=nanmean(a(:,[5,6,end]),2);  %for checking      oo=a(:,[5,6,end]);
        thres_vis=[a(:,1:4),round([ann,spring,summer,autumn,winter],0)];
        
        xname=sort_nat_name2{ix};
        headername=xname(1:length(xname)-20);
        tailname=xname(length(xname)-17:end);
        midname='年四季';
        a_newname=strcat(headername,midname,tailname);
        
        str1={'tag','lon','lat','ht','年','春','夏','秋','冬'};
        fname2=a_newname;
        fid2=fopen([subpath,fname2],'w+');%w+是重写内容
        for ii=1:size(str1,2)
            fprintf(fid2,'%s\t',str1{1,ii});%'%s\t'表明输出字符串后加一个tab键，好跟后面的字符串分开
        end
        fprintf(fid2,'\n');
        fclose(fid2);
        dlmwrite([subpath,fname2],thres_vis,'delimiter','\t','-append');
%-append表示末尾追加，在第一行str1换行后追加thres_vis（thres_vis已经是所有台站结果了）   
    end  
end

%%  5.计算不同条件下四季或者逐月的medianVIS均值--区域阈值 
%{
%基于区域数据：91台站*12个月。计算冬季，某站12月，1月，2月，算数平均作为该站冬季的阈值，
%其他季节亦如此，如此就有91台站*四季；再以此为依据，某季节的91个阈值做算术平均作为该季节区域阈值

程序思路：某季节的91个阈值做算术平均作为该季节或者月份的区域阈值
%}
clear,clc
list=dir('Z:\data\obs\2001-2021pre-t-rh\pro\matfile\2023*');
%手动点击list看看，有name对应三个文件夹，folder为上级路径{'Z:\data\obs\2001-2021pre-t-rh\pro\matfile'}
nf=size(list,1); %获取文件夹的个数  
matpath=list(1).folder;  

%-----------------
fname1='基于中位数-2O23文件排序-年季节2016-2020的VIS区域阈值.txt';
fid1=fopen([matpath,'\',fname1],'w+');
fprintf(fid1,'%s\r\n','第一行是某降水等级下所有台站一起在某年或季节VIS的区域均值，第二行是max');%'%s\t'表明输出字符串后加一个tab键，好跟后面的字符串分开
fprintf(fid1,'%s\r\n','年  春季   夏季   秋季  冬季');
fclose(fid1);
fname2='基于中位数-2O23文件排序-monthly2016-2020的VIS区域阈值.txt';
fid2=fopen([matpath,'\',fname2],'w+');
fprintf(fid2,'%s\r\n','第一行是某降水等级下所有台站一起在某月份medianVIS的区域均值，第二行是max');%'%s\t'表明输出字符串后加一个tab键，好跟后面的字符串分开
fprintf(fid2,'%s\r\n','1月    2月  3月  4月  5月  6月  7月  8月  9月  10月 11月 12月');
fclose(fid2);
%-----------------

for i=1:nf
    %subfoldername=list.name{i,:};%wrong wrong wrong
    subfoldername=list(i).name;
    subpath=strcat(matpath,'\',subfoldername,'\');
    sublist1=dir([subpath,'中位数*四季*.txt']); 
    sublist2=dir([subpath,'中位数*逐月*.txt']); 
    sort_nat_name1=sort_nat({sublist1.name}); %使用sort_nat函数对文件排序
    sort_nat_name2=sort_nat({sublist2.name}); %
    %=======
    %   seasonal
    %=======
    %{
    formats1=['%d',repmat('%f',[1,3]),repmat('%s',[1,5])];
    for ix=1:size(sublist1,1)  %txt文件的个数 
        fid1=fopen([subpath,sublist1(ix).name]);
        a=textscan(fid1,formats1,'headerlines',1,'delimiter','\t','CollectOutput',1);%        
        fclose(fid1);
        dat3=a{1,3};%91*5cell
        aimdat=cell2mat(dat3);%总出错，不用repmat('%s',[1,5]，文件里的NaN被自动转换为0
        
    end
    %}
    for ix=1:size(sublist1,1)  %txt文件的个数 
        a=dlmread([subpath,sort_nat_name1{ix}],'\t',1,0);%用排序后的文件名读入
        %a=dlmread([subpath,sublist1(ix).name],'\t',1,0);%左上角的起始坐标是R=0,C=0 ，得到91*9的double数据  
        %a2=importdata([subpath,sublist1(ix).name],'\t',1);%1表示取第一行以下，得到一个结构数据
        %a is tag-lat-lon-ht-年四季
        season_mins=nanmean(a(:,5:end));%nanmin(a(:,5:end));
        season_maxs=nanmax(a(:,5:end));
        aout=round([season_mins;season_maxs],0);
        fid1=fopen([matpath,'\',fname1],'a');%以追加模式打开  
        fprintf(fid1,'%s\r\n',sort_nat_name1{ix});%'%s\t'表明输出字符串后加一个tab键，好跟后面的字符串分开
        fprintf(fid1,'%d  %d  %d  %d  %d\r\n',aout');
        fclose(fid1);
        disp([i,ix])  %确定它遍历了三个文件夹
%         dlmwrite([matpath,'\',fname1],[sublist1(ix).name,'\n'],'-append');%dlmwrite自动分割字符串          
%         dlmwrite([matpath,'\',fname1],aout,'delimiter','\t','-append');        
    end   
    
    %=======
    %   monthly
    %=======    
    for iy=1:size(sublist1,1)  %txt文件的个数 
        b=dlmread([subpath,sort_nat_name2{iy}],'\t',1,0);%用排序后的文件名读入
        %b=dlmread([subpath,sublist2(iy).name],'\t',1,0);%左上角的起始坐标是R=0,C=0 ，得到91*9的double数据  
        %b is tag-lat-lon-ht-1月~12月
        mon_mins=nanmean(b(:,5:end));%
        mon_maxs=nanmax(b(:,5:end));
        bout=round([mon_mins;mon_maxs],0);
        fid2=fopen([matpath,'\',fname2],'a');%以追加模式打开  
        fprintf(fid2,'%s\r\n',sort_nat_name2{iy});%'%s\t'表明输出字符串后加一个tab键，好跟后面的字符串分开
        fprintf(fid2,[repmat('%d  ',[1,11]),'%d\r\n'],bout');
        fclose(fid2);
        disp([i,iy])  %确定它遍历了三个文件夹      
    end     
    
end

%%  6.基于中位数，去中文名，转写数据：四季或者逐月
%只需针对年四季的数据转写即可，逐月数据的结果不变
clear,clc
%转写为第一行台站，第二行纬度，第三行经度，第四行到第尾行依次是年、春夏秋冬

list=dir('Z:\data\obs\2001-2021pre-t-rh\pro\matfile\2023*');
%手动点击list看看，有name对应三个文件夹，folder为上级路径{'Z:\data\obs\2001-2021pre-t-rh\pro\matfile'}
nf=size(list,1); %获取文件夹的个数  
matpath=list(1).folder;  
outpath='Z:\data\obs\2001-2021pre-t-rh\pro\matfile\forGrads\';
%-----------------

%-----------------

for i=1:nf
    %subfoldername=list.name{i,:};%wrong wrong wrong
    subfoldername=list(i).name;
    subpath=strcat(matpath,'\',subfoldername,'\');
    sublist1=dir([subpath,'中位数from*四季*.txt']); 
    sublist2=dir([subpath,'中位数from*逐月*.txt']); 
    %=======
    %   seasonal
    %=======

    for ix=1:size(sublist1,1)  %txt文件的个数 
        a=dlmread([subpath,sublist1(ix).name],'\t',1,0);%左上角的起始坐标是R=0,C=0 ，得到91*9的double数据  
        %a2=importdata([subpath,sublist1(ix).name],'\t',1);%1表示取第一行以下，得到一个结构数据
        %a is tag-lon-lat-ht-年四季
        %sublist1(ix).name is 'from0.1to1.5mm各站年四季vis阈值2016-2020.txt'
        aout=[a(:,1),a(:,3),a(:,2),a(:,5:end)];%aout is tag-lat-lon-年四季 
        aout=aout'; %转写为第一行台站，第二行纬度，第三行经度，第四行到第尾行依次是年、春夏秋冬
        aout(isnan(aout))=-999;%32766; %0406修改添加
        
        %按这种fmt1输出91个台站后不换行哦
        %fmt1=repmat('%-8d',[1 91]);%默认右对齐，加个-号变为左对齐
        %fmt2=repmat('%-8.2f',[1 91]);  
        fmt1=[repmat('%-8d',[1 90]),'%-8d\r\n'];%一个fmt1末尾换行
        fmt2=[repmat('%-8.2f',[1 90]),'%-8.2f\r\n'];%一个fmt2末尾换行   
        
        %grads不认中文
        xname=sublist1(ix).name;
        headername=xname(4:length(xname)-23);%去掉质控后三个字
        tailname=xname(length(xname)-12:end);
        midname='_ANNseasonsVISthres';
        a_newname=strcat(headername,midname,tailname);
        
        fid1=fopen([outpath,a_newname],'w+');%以覆盖模式打开
        %fid1=fopen([outpath,'trans_',sublist1(ix).name],'w+');%以覆盖模式打开
        %fid1=fopen([subpath,'trans_',sublist1(ix).name],'w+');%以覆盖模式打开,NaN不做替换  
        %fprintf(fid1,[fmt1,fmt2,fmt2,fmt1,fmt1,fmt1,fmt1,'%8d\r\n'],aout');%old
        fprintf(fid1,[fmt1,fmt2,fmt2,fmt1,fmt1,fmt1,fmt1,fmt1],aout');
        fclose(fid1);
        disp([i,ix])  %确定它遍历了三个文件夹    
    end    
    
    
    %=======
    %   monthly
    %=======    
    for iy=1:size(sublist1,1)  %txt文件的个数 
        b=dlmread([subpath,sublist2(iy).name],'\t',1,0);%左上角的起始坐标是R=0,C=0 ，得到91*9的double数据  
        %b is tag-lat-lon-ht-1月~12月
        %bout=[b(:,1),b(:,3),b(:,2),b(:,5:end)];%bout is tag-lat-lon-1到12月
        bout=[b(:,1),b(:,3),b(:,2),round(b(:,5:end),0)];%bout is tag-lat-lon-1到12月，必须有round(b(:,5:end),0)
        bout=bout'; %转写为第一行台站，第二行纬度，第三行经度，第四行到第尾行依次是年、春夏秋冬
        bout(isnan(bout))=-999;%32766; %0406修改添加
        
       %grads不认中文
        xname=sublist1(iy).name; %attention sublist1(iy).name
        headername=xname(4:length(xname)-23);%去掉质控后三个字
        tailname=xname(length(xname)-12:end);
        midname='_MonthlyVISthres';
        b_newname=strcat(headername,midname,tailname);
        
        fid2=fopen([outpath,b_newname],'w+');%以覆盖模式打开
        %fid2=fopen([outpath,'trans_',sublist2(iy).name],'w+');%以覆盖模式打开
        %fid2=fopen([subpath,'trans_',sublist2(iy).name],'w+');%以覆盖模式打开,NaN不做替换
        %fprintf(fid2,[fmt1,fmt2,fmt2,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,'%8d\r\n'],bout');%old
        fprintf(fid2,[fmt1,fmt2,fmt2,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1],bout');
        fclose(fid2);
        disp([i,iy])  %确定它遍历了三个文件夹      
    end     
    
end

%%  ANNorMonthly手动更换   英文文件名，批量写ctl
%{
%需要了解的信息，自己查看文件的名称等
pathin_obsgrd='E:\fortran\stn2grd\';%obs数据转二进制后存放的路径
pathin_gridgrd='E:\fortran\stn2grd\netgrids\';%网格点转二进制后存放的路径
header_gridgrd='grids_';%网格点转二进制后的文件打头
tail_grd='.txt.grd';%转二进制时直接在txt后面加了.grd

%}

clear,clc
%----待写ctl的文件名----
pathin='Z:\data\obs\2001-2021pre-t-rh\pro\matfile\forGrads\';%其中的seasonal.txt和monthly.txt
%手动修改以下两处
filelist='monthly.txt'; %monthly.txt  or  seasonal.txt
keywords='monthly';% seasonal   or monthly

fid=fopen([pathin,filelist],'r');
x=textscan(fid,'%s');%x is 1*1 cell,
xx=x{1,1}; %x{1,1}才是13*1cell,包括了13个文件名
fclose(fid); 

%----导入ctl样本------
ctlpathin='Z:\STUDY-Storage-box\grads\GDScripts\medianVIS\';
pathout=ctlpathin;

%查看aa后，确定站点数据的ctl，需要修改第二、六行
fid0=fopen([ctlpathin,keywords,'VIS.ctl'],'r');%
a=textscan(fid0,'%s');%x is 1*1 cell,
aa=a{1,1}; %a{1,1}才是22*1cell, 
fclose(fid0); 

%查看cc后，确定网格点的ctl，需要修改第二行
fid0=fopen([ctlpathin,keywords,'VISgrid.ctl'],'r');
c=textscan(fid0,'%s');%x is 1*1 cell,
cc=c{1,1}; %c{1,1}才是36*1cell,
fclose(fid0); 

% % %========================================
% % %根据待写ctl的文件名xx进行对站点数据和网格点ctl的修改
% % %========================================
%需要了解的信息，自己查看文件的名称等，这个根据存放路径而修改
pathin_obsgrd='Z:\STUDY-Storage-box\grads\GDScripts\medianVIS\grd\';%obs数据转二进制后存放的路径
pathin_gridgrd='Z:\STUDY-Storage-box\grads\GDScripts\medianVIS\grd\';%网格点转二进制后存放的路径
header_obsgrd='';%obs数据转二进制后的文件名打头
header_gridgrd='grids_';%网格点转二进制后的文件打头
tail_grd='.txt.grd';%转二进制时直接在txt后面加了.grd

for i=1:size(xx,1)
    fname=xx{i};%eg: from0.1to0.4mm_ANNseasonsVISthres2016-2020.txt
    newname=fname(1:length(fname)-4);%去尾后：from0.1to0.4mm_ANNseasonsVISthres2016-2020
    tmpaa=aa;
    tmpaa{2}=strcat(pathin_obsgrd,header_obsgrd,newname,tail_grd); %用于写入ctl中，识别台站数据的二进制
    tmpaa{6}=strcat(pathout,newname,'.map');
    %%%对照样本的ctl组合一下读入的ctl数据
    newaa{1}=[tmpaa{1},'  ',tmpaa{2}];
    newaa{2}=[tmpaa{3},'  ',tmpaa{4}];
    newaa{3}=[tmpaa{5},'  ',tmpaa{6}];
    newaa{4}=[tmpaa{7},'  ',tmpaa{8}];
    newaa{5}=[tmpaa{9},'  ',tmpaa{10}];
    newaa{6}=[tmpaa{11},'  ',tmpaa{12},'  ',tmpaa{13},'  ',tmpaa{14},'  ',tmpaa{15}];
    newaa{7}=[tmpaa{16},'  ',tmpaa{17}];
    newaa{8}=[tmpaa{18},'  ',tmpaa{19},'  ',tmpaa{20},'  ',tmpaa{21}];
    newaa{9}=tmpaa{22};    
    %输出obs二进制对应的ctl
    fid2=fopen([pathout,newname,'.ctl'],'w+');
    fprintf(fid2,'%s\r\n',newaa{:});
    fclose(fid2); 
 
    tmpcc=cc;
    tmpcc{2}=strcat(pathin_gridgrd,header_gridgrd,newname,tail_grd);
    newcc{1}=[tmpcc{1},'  ',tmpcc{2}];
    newcc{2}=[tmpcc{3},'  ',tmpcc{4}];
    newcc{3}=[tmpcc{5},'  ',tmpcc{6},'  ',tmpcc{7},'  ',tmpcc{8}];
    newcc{4}=[tmpcc{9},'  ',tmpcc{10},'  ',tmpcc{11},'  ',tmpcc{12},'  ',tmpcc{13}];
    newcc{5}=[tmpcc{14},'  ',tmpcc{15},'  ',tmpcc{16},'  ',tmpcc{17},'  ',tmpcc{18}];
    newcc{6}=[tmpcc{19},'  ',tmpcc{20},'  ',tmpcc{21},'  ',tmpcc{22},'  ',tmpcc{23}];
    newcc{7}=[tmpcc{24},'  ',tmpcc{25},'  ',tmpcc{26},'  ',tmpcc{27},'  ',tmpcc{28}];
    newcc{8}=[tmpcc{29},'  ',tmpcc{30}];
    newcc{9}=[tmpcc{31},'  ',tmpcc{32},'  ',tmpcc{33},'  ',tmpcc{34},'  ',tmpcc{35}];
    newcc{10}=[tmpcc{36}]; 
    %输出grid二进制对应的ctl
    fid3=fopen([pathout,newname,'grid.ctl'],'w+');
    fprintf(fid3,'%s\r\n',newcc{:});
    fclose(fid3); 
    
end