%{
llyu@04.24   
��˼·�����㲻ͬ��ˮ�ȼ��¸����ܼ���������λ���γɵ�91վ*12����
�Դ�Ϊ�������ֱ������ͼ��ڣ��õ�91վ*5���ļ�������(��)ֵ=��Ӧ�·���ӳ���3(12)
������ֵ�ķ�����nanmean(91վ*12����)  ��   nanmean(91վ*5���ļ�)
%}
%% 1.cal-hourly��ˮ�ֳ�С�д��¸�վmedianVIS--����
%VIS�Զ��۲�ʼ��2016
%%tag-yyyy-mm-dd-hh-t(���϶�)-pre��mm��-rh(%)-vis(m)
%��ϸ�۲����ݣ�ʱ�̴�21ʱ��20ʱ�����˽�ˮΪ��ʱ���ݣ���������һ���Ĵ�
%�۲���ĸ�ʱ��Ϊ02/08/14/20����ȻΪBJT
%˼·���ȵȼ���������ȡ��λ��
clear,clc
pathin='Z:\data\obs\2001-2021pre-t-rh\pro\ʱ�������\';
pathout='Z:\data\obs\2001-2021pre-t-rh\pro\';
matpath='Z:\data\obs\2001-2021pre-t-rh\pro\matfile\';
pltpath='Z:\data\obs\2001-2021pre-t-rh\fig\';
load([matpath,'stn91_classify.mat']);
%��������station_info��station_types�� 'plain_stns','hill_stns','mountain_stns'
xstation_types={'ƽԭվ','����վ','��ɽվ'};%��Ӧ{'plain','hill','mountain'}

PREclass.str{1,1}='��ʱ��ˮ����ΪС�ꡢ���ꡢ����ͱ����ĸ��ȼ�';
PREclass.str{2,1}='��ˮ�ȼ���׼-����GB/T25892-2012';
PREclass.cn={'����','С��','����','����','����'};
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
        %%tag-yyyy-mm-dd-hh-t(���϶�)��6��-pre��mm����7��-rh(%)��8��-vis(m)��9��
        %VIS�Զ��۲�ʼ��2016,��ȡ1��1��0ʱ��  ����ʱ������
        line_2016=find(raw(:,2)==2016 & raw(:,3)==1 & raw(:,4)==1 & raw(:,5)==0); 
        line_end=find(raw(:,2)==2020 & raw(:,3)==12 & raw(:,4)==31 & raw(:,5)==23);
        sample2x=raw(line_2016:line_end,:);%visû��ȡ2021�꣬��Ϊ2021���������ڼ�����        
        %sample2=sample2x(sample2x(:,7)>=0.1,:);%�о�vis�����ּ���ˮ��ʱ������,�������ֿ��Բ�����

        vec_vis=zeros(1,12)*nan;%minvis[1to12]          
        for im=1:12  %�·�
            %��ȡ�ȼ�����---���ݽ�ˮ��
            class_sample=sample2x((sample2x(:,7)>=PREclass.below(ic) & sample2x(:,7)<=PREclass.up(ic)),:);
            %tag-yyyy-mm-dd-hh-t(���϶�)��6��-pre��mm����7��-rh(%)��8��-vis(m)��9��            
            
            %{
            %��ȥ�쳣ֵ---������׼��
            xvec=class_sample(:,9); %�����������
            xave=nanmean(xvec);   xstd=nanstd(xvec); %��ֵ�;�����
            up_val=xave+3*xstd;  below_val=xave-3*xstd; %���޺�����
            inx=find(xvec<=up_val & xvec>=below_val); %����ֵ��Ӧ������
            tmpx=class_sample(inx,:); % ����ֵ
            iny=find(xvec>up_val | xvec<below_val); %�쳣ֵ��Ӧ������
            tmpyy=class_sample(iny,:); % to check �쳣ֵ  
            %���ȡ����
            sample_for_vis=tmpx(tmpx(:,3)==im,:);%�о�vis�����н�ˮ��ʱ������            
            %}            
            %���ȡ����
            sample_for_vis=class_sample(class_sample(:,3)==im,:);%�о�vis�����н�ˮ��ʱ������
            %==============================================
            %-----------��ĳ�ȼ���ˮĳ�¶�Ӧ��vis��Сֵ-----------            
            if ~isempty(sample_for_vis)
                %minvis=nanmin(sample_for_vis(:,9));%���ֽ�ˮ��Ӧ��vis��ֵ
                minvis=nanmedian(sample_for_vis(:,9));%���ֽ�ˮ��Ӧ��vis��ֵ
            else
                minvis=nan;
            end
            vec_vis(1,im)=round(minvis,0);   %vec_vis(1,im)=minvis;

        end %��������ֵ im
        xthres_vis(id,:)=[tagsinfo(id,:),vec_vis];
        %xthres_visһ��16�У�tag-lon-lat-ht-minvis[1-12��]
    end %ĳ��̨վѭ��  id
    eval(['minvis_',num2str(itype),'=xthres_vis;']);%tag-lon-lat-ht-minvis[1-12��]
end %̨վ���� itype
thres_vis=[minvis_1;minvis_2;minvis_3];%����̨վ�Ľ����ֱ���ӣ�û����׷����ʽ
str1={'tag','lon','lat','ht','1��','2��','3��','4��','5��','6��','7��','8��','9��','10��','11��','12��'};

%str1={'tag','lon','lat','ht','��','��','��','��','��'};
%str1={'tag','lon','lat','ht','��','��','��','��'};
%dlmwrite([matpath,'��վ�ļ�vis��ֵ.txt'],str1{:},'delimiter','\t');
%'precision','%.4f',

fname2=strcat('��λ��from',num2str(PREclass.below(ic)),'to',num2str(PREclass.up(ic)),'mm��վ����vis��ֵ2016-2020.txt');
fid2=fopen([matpath,fname2],'w+');%w+����д����
for ii=1:size(str1,2)
fprintf(fid2,'%s\t',str1{1,ii});%'%s\t'��������ַ������һ��tab�����ø�������ַ����ֿ�
end
fprintf(fid2,'\n');
fclose(fid2);
dlmwrite([matpath,fname2],thres_vis,'delimiter','\t','-append');
%-append��ʾĩβ׷�ӣ��ڵ�һ��str1���к�׷��thres_vis��thres_vis�Ѿ�������̨վ����ˣ�
%������������������������������������������������
end %��ˮ�ȼ� ic

%% 2.cal-hourly��ˮ�ֳ�6���ȼ���վmedianVIS--����
%VIS�Զ��۲�ʼ��2016
%%tag-yyyy-mm-dd-hh-t(���϶�)-pre��mm��-rh(%)-vis(m)
%��ϸ�۲����ݣ�ʱ�̴�21ʱ��20ʱ�����˽�ˮΪ��ʱ���ݣ���������һ���Ĵ�
%�۲���ĸ�ʱ��Ϊ02/08/14/20����ȻΪBJT
%���ֽ�ˮʱ�̶�Ӧ��VIS��ʷ���ֵ
%˼·���ȵȼ���������ȡ��λ��

clear,clc
pathin='Z:\data\obs\2001-2021pre-t-rh\pro\ʱ�������\';
pathout='Z:\data\obs\2001-2021pre-t-rh\pro\';
matpath='Z:\data\obs\2001-2021pre-t-rh\pro\matfile\';
pltpath='Z:\data\obs\2001-2021pre-t-rh\fig\';
load([matpath,'stn91_classify.mat']);
%��������station_info��station_types�� 'plain_stns','hill_stns','mountain_stns'
xstation_types={'ƽԭվ','����վ','��ɽվ'};%��Ӧ{'plain','hill','mountain'}

PREclass.str{1,1}='��ʱ��ˮ����Ϊ�����ȼ�';
PREclass.str{2,1}='��ˮ�ȼ���׼-�Զ���';
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
        %%tag-yyyy-mm-dd-hh-t(���϶�)��6��-pre��mm����7��-rh(%)��8��-vis(m)��9��
        %VIS�Զ��۲�ʼ��2016,��ȡ1��1��0ʱ��  ����ʱ������
        line_2016=find(raw(:,2)==2016 & raw(:,3)==1 & raw(:,4)==1 & raw(:,5)==0); 
        line_end=find(raw(:,2)==2020 & raw(:,3)==12 & raw(:,4)==31 & raw(:,5)==23);
        sample2x=raw(line_2016:line_end,:);%visû��ȡ2021�꣬��Ϊ2021���������ڼ�����        
        %sample2=sample2x(sample2x(:,7)>=0.1,:);%�о�vis�����ּ���ˮ��ʱ������,�������ֿ��Բ�����

        vec_vis=zeros(1,12)*nan;%minvis[1to12]          
        for im=1:12  %�·�
            %��ȡ�ȼ�����---���ݽ�ˮ��
            class_sample=sample2x((sample2x(:,7)>=PREclass.below(ic) & sample2x(:,7)<=PREclass.up(ic)),:);
            %tag-yyyy-mm-dd-hh-t(���϶�)��6��-pre��mm����7��-rh(%)��8��-vis(m)��9��            
            
            %{
            %��ȥ�쳣ֵ---������׼��
            xvec=class_sample(:,9); %�����������
            xave=nanmean(xvec);   xstd=nanstd(xvec); %��ֵ�;�����
            up_val=xave+3*xstd;  below_val=xave-3*xstd; %���޺�����
            inx=find(xvec<=up_val & xvec>=below_val); %����ֵ��Ӧ������
            tmpx=class_sample(inx,:); % ����ֵ
            iny=find(xvec>up_val | xvec<below_val); %�쳣ֵ��Ӧ������
            tmpyy=class_sample(iny,:); % to check �쳣ֵ             
            %���ȡ����
            sample_for_vis=tmpx(tmpx(:,3)==im,:);%�о�vis�����н�ˮ��ʱ������
            %}
            
            sample_for_vis=class_sample(class_sample(:,3)==im,:);%�о�vis�����н�ˮ��ʱ������
            %==============================================
            %-----------��ĳ�ȼ���ˮĳ�¶�Ӧ��vis��Сֵ-----------            
            if ~isempty(sample_for_vis)
                minvis=nanmedian(sample_for_vis(:,9));%���ֽ�ˮ��Ӧ��vis��ֵ
                %minvis=nanmin(sample_for_vis(:,9));%���ֽ�ˮ��Ӧ��vis��ֵ
            else
                minvis=nan;
            end
            vec_vis(1,im)=round(minvis,0);   %vec_vis(1,im)=minvis;

        end %��������ֵ im
        xthres_vis(id,:)=[tagsinfo(id,:),vec_vis];
        %xthres_visһ��16�У�tag-lon-lat-ht-minvis[1-12��]
    end %ĳ��̨վѭ��  id
    eval(['minvis_',num2str(itype),'=xthres_vis;']);%tag-lon-lat-ht-minvis[1-12��]
end %̨վ���� itype
thres_vis=[minvis_1;minvis_2;minvis_3];%����̨վ�Ľ����ֱ���ӣ�û����׷����ʽ
str1={'tag','lon','lat','ht','1��','2��','3��','4��','5��','6��','7��','8��','9��','10��','11��','12��'};

%str1={'tag','lon','lat','ht','��','��','��','��','��'};
%str1={'tag','lon','lat','ht','��','��','��','��'};
%dlmwrite([matpath,'��վ�ļ�vis��ֵ.txt'],str1{:},'delimiter','\t');
%'precision','%.4f',

fname2=strcat('��λ��from',num2str(PREclass.below(ic)),'to',num2str(PREclass.up(ic)),'mm��վ����vis��ֵ2016-2020.txt');
fid2=fopen([matpath,fname2],'w+');%w+����д����
for ii=1:size(str1,2)
fprintf(fid2,'%s\t',str1{1,ii});%'%s\t'��������ַ������һ��tab�����ø�������ַ����ֿ�
end
fprintf(fid2,'\n');
fclose(fid2);
dlmwrite([matpath,fname2],thres_vis,'delimiter','\t','-append');
%-append��ʾĩβ׷�ӣ��ڵ�һ��str1���к�׷��thres_vis��thres_vis�Ѿ�������̨վ����ˣ�
%������������������������������������������������
end %��ˮ�ȼ� ic

%% 3.cal-hourly��ˮ���ֵȼ���վmedianVIS--����
%VIS�Զ��۲�ʼ��2016
%%tag-yyyy-mm-dd-hh-t(���϶�)-pre��mm��-rh(%)-vis(m)
%��ϸ�۲����ݣ�ʱ�̴�21ʱ��20ʱ�����˽�ˮΪ��ʱ���ݣ���������һ���Ĵ�
%�۲���ĸ�ʱ��Ϊ02/08/14/20����ȻΪBJT
%���ֽ�ˮʱ�̶�Ӧ��VIS��ʷ���ֵ
%˼·���ȵȼ���������ȡ��λ��
clear,clc
pathin='Z:\data\obs\2001-2021pre-t-rh\pro\ʱ�������\';
pathout='Z:\data\obs\2001-2021pre-t-rh\pro\';
matpath='Z:\data\obs\2001-2021pre-t-rh\pro\matfile\';
pltpath='Z:\data\obs\2001-2021pre-t-rh\fig\';
load([matpath,'stn91_classify.mat']);
%��������station_info��station_types�� 'plain_stns','hill_stns','mountain_stns'
xstation_types={'ƽԭվ','����վ','��ɽվ'};%��Ӧ{'plain','hill','mountain'}

PREclass.str{1,1}='���ֵȼ�';
PREclass.str{2,1}='���ڵ���0.1��0.5mm';
PREclass.cn={'���ڵ���0.1','���ڵ���0.5'};
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
        %%tag-yyyy-mm-dd-hh-t(���϶�)��6��-pre��mm����7��-rh(%)��8��-vis(m)��9��
        %VIS�Զ��۲�ʼ��2016,��ȡ1��1��0ʱ��  ����ʱ������
        line_2016=find(raw(:,2)==2016 & raw(:,3)==1 & raw(:,4)==1 & raw(:,5)==0); 
        line_end=find(raw(:,2)==2020 & raw(:,3)==12 & raw(:,4)==31 & raw(:,5)==23);
        sample2x=raw(line_2016:line_end,:);%visû��ȡ2021�꣬��Ϊ2021���������ڼ�����        
        %sample2=sample2x(sample2x(:,7)>=0.1,:);%�о�vis�����ּ���ˮ��ʱ������,�������ֿ��Բ�����

        vec_vis=zeros(1,12)*nan;%minvis[1to12]          
        for im=1:12  %�·�
            %��ȡ�ȼ�����---���ݽ�ˮ��
            class_sample=sample2x((sample2x(:,7)>=PREclass.below(ic) & sample2x(:,7)<=PREclass.up(ic)),:);
            %tag-yyyy-mm-dd-hh-t(���϶�)��6��-pre��mm����7��-rh(%)��8��-vis(m)��9��            
            
            %{
            %��ȥ�쳣ֵ---������׼��
            xvec=class_sample(:,9); %�����������
            xave=nanmean(xvec);   xstd=nanstd(xvec); %��ֵ�;�����
            up_val=xave+3*xstd;  below_val=xave-3*xstd; %���޺�����
            inx=find(xvec<=up_val & xvec>=below_val); %����ֵ��Ӧ������
            tmpx=class_sample(inx,:); % ����ֵ
            iny=find(xvec>up_val | xvec<below_val); %�쳣ֵ��Ӧ������
            tmpyy=class_sample(iny,:); % to check �쳣ֵ              
            %���ȡ����
            sample_for_vis=tmpx(tmpx(:,3)==im,:);%�о�vis�����н�ˮ��ʱ������
            %}
            
            sample_for_vis=class_sample(class_sample(:,3)==im,:);%�о�vis�����н�ˮ��ʱ������
            %==============================================
            %-----------��ĳ�ȼ���ˮĳ�¶�Ӧ��vis��Сֵ-----------                
            if ~isempty(sample_for_vis)
                minvis=nanmedian(sample_for_vis(:,9));%���ֽ�ˮ��Ӧ��vis��ֵ
                %minvis=nanmin(sample_for_vis(:,9));%���ֽ�ˮ��Ӧ��vis��ֵ
            else
                minvis=nan;
            end
            vec_vis(1,im)=round(minvis,0);   %vec_vis(1,im)=minvis;

        end %��������ֵ im
        xthres_vis(id,:)=[tagsinfo(id,:),vec_vis];
        %xthres_visһ��16�У�tag-lon-lat-ht-minvis[1-12��]
    end %ĳ��̨վѭ��  id
    eval(['minvis_',num2str(itype),'=xthres_vis;']);%tag-lon-lat-ht-minvis[1-12��]
end %̨վ���� itype
thres_vis=[minvis_1;minvis_2;minvis_3]; %����̨վ�Ľ����ֱ���ӣ�û����׷����ʽ
str1={'tag','lon','lat','ht','1��','2��','3��','4��','5��','6��','7��','8��','9��','10��','11��','12��'};

%str1={'tag','lon','lat','ht','��','��','��','��','��'};
%str1={'tag','lon','lat','ht','��','��','��','��'};
%dlmwrite([matpath,'��վ�ļ�vis��ֵ.txt'],str1{:},'delimiter','\t');
%'precision','%.4f',

fname2=strcat('��λ��from',num2str(PREclass.below(ic)),'to',num2str(PREclass.up(ic)),'mm��վ����vis��ֵ2016-2020.txt');
fid2=fopen([matpath,fname2],'w+');%w+����д����
for ii=1:size(str1,2)
fprintf(fid2,'%s\t',str1{1,ii});%'%s\t'��������ַ������һ��tab�����ø�������ַ����ֿ�
end
fprintf(fid2,'\n');
fclose(fid2);
dlmwrite([matpath,fname2],thres_vis,'delimiter','\t','-append');
%-append��ʾĩβ׷�ӣ��ڵ�һ��str1���к�׷��thres_vis��thres_vis�Ѿ�������̨վ����ˣ�
%������������������������������������������������
end %��ˮ�ȼ� ic  

%% 4.���㲻ͬ�����¸�վ����ļ�����ֵ������91̨վ*12����medianVIS--���ļ� 
%�����������ݣ�91̨վ*12���¡����㶬����ĳվ12�£�1�£�2�£�����ƽ����Ϊ��վ��������ֵ��
%������������ˣ���˾���91̨վ*�ļ������Դ�Ϊ���ݣ�ĳ���ڵ�91����ֵ������ƽ����Ϊ�ü���������ֵ
clear,clc
list=dir('Z:\data\obs\2001-2021pre-t-rh\pro\matfile\2023*');
%�ֶ����list��������name��Ӧ�����ļ��У�folderΪ�ϼ�·��{'Z:\data\obs\2001-2021pre-t-rh\pro\matfile'}
nf=size(list,1); %��ȡ�ļ��еĸ���  
matpath=list(1).folder;  

for i=1:nf
    %subfoldername=list.name{i,:};%wrong wrong wrong
    subfoldername=list(i).name;
    subpath=strcat(matpath,'\',subfoldername,'\');
     sublist2=dir([subpath,'��λ��*����*.txt']); 
    sort_nat_name2=sort_nat({sublist2.name}); %ʹ��sort_nat�������ļ�����
    %=======
    %   seasonal
    %=======
    for ix=1:size(sublist2,1)  %txt�ļ��ĸ��� 
        a=dlmread([subpath,sort_nat_name2{ix}],'\t',1,0);%���������ļ�������
        %a=dlmread([subpath,sublist1(ix).name],'\t',1,0);%���Ͻǵ���ʼ������R=0,C=0 ���õ�91*9��double����  
        %a2=importdata([subpath,sublist1(ix).name],'\t',1);%1��ʾȡ��һ�����£��õ�һ���ṹ����
        %a is tag-lat-lon-ht-1��12����
        
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
        midname='���ļ�';
        a_newname=strcat(headername,midname,tailname);
        
        str1={'tag','lon','lat','ht','��','��','��','��','��'};
        fname2=a_newname;
        fid2=fopen([subpath,fname2],'w+');%w+����д����
        for ii=1:size(str1,2)
            fprintf(fid2,'%s\t',str1{1,ii});%'%s\t'��������ַ������һ��tab�����ø�������ַ����ֿ�
        end
        fprintf(fid2,'\n');
        fclose(fid2);
        dlmwrite([subpath,fname2],thres_vis,'delimiter','\t','-append');
%-append��ʾĩβ׷�ӣ��ڵ�һ��str1���к�׷��thres_vis��thres_vis�Ѿ�������̨վ����ˣ�   
    end  
end

%%  5.���㲻ͬ�������ļ��������µ�medianVIS��ֵ--������ֵ 
%{
%�����������ݣ�91̨վ*12���¡����㶬����ĳվ12�£�1�£�2�£�����ƽ����Ϊ��վ��������ֵ��
%������������ˣ���˾���91̨վ*�ļ������Դ�Ϊ���ݣ�ĳ���ڵ�91����ֵ������ƽ����Ϊ�ü���������ֵ

����˼·��ĳ���ڵ�91����ֵ������ƽ����Ϊ�ü��ڻ����·ݵ�������ֵ
%}
clear,clc
list=dir('Z:\data\obs\2001-2021pre-t-rh\pro\matfile\2023*');
%�ֶ����list��������name��Ӧ�����ļ��У�folderΪ�ϼ�·��{'Z:\data\obs\2001-2021pre-t-rh\pro\matfile'}
nf=size(list,1); %��ȡ�ļ��еĸ���  
matpath=list(1).folder;  

%-----------------
fname1='������λ��-2O23�ļ�����-�꼾��2016-2020��VIS������ֵ.txt';
fid1=fopen([matpath,'\',fname1],'w+');
fprintf(fid1,'%s\r\n','��һ����ĳ��ˮ�ȼ�������̨վһ����ĳ��򼾽�VIS�������ֵ���ڶ�����max');%'%s\t'��������ַ������һ��tab�����ø�������ַ����ֿ�
fprintf(fid1,'%s\r\n','��  ����   �ļ�   �＾  ����');
fclose(fid1);
fname2='������λ��-2O23�ļ�����-monthly2016-2020��VIS������ֵ.txt';
fid2=fopen([matpath,'\',fname2],'w+');
fprintf(fid2,'%s\r\n','��һ����ĳ��ˮ�ȼ�������̨վһ����ĳ�·�medianVIS�������ֵ���ڶ�����max');%'%s\t'��������ַ������һ��tab�����ø�������ַ����ֿ�
fprintf(fid2,'%s\r\n','1��    2��  3��  4��  5��  6��  7��  8��  9��  10�� 11�� 12��');
fclose(fid2);
%-----------------

for i=1:nf
    %subfoldername=list.name{i,:};%wrong wrong wrong
    subfoldername=list(i).name;
    subpath=strcat(matpath,'\',subfoldername,'\');
    sublist1=dir([subpath,'��λ��*�ļ�*.txt']); 
    sublist2=dir([subpath,'��λ��*����*.txt']); 
    sort_nat_name1=sort_nat({sublist1.name}); %ʹ��sort_nat�������ļ�����
    sort_nat_name2=sort_nat({sublist2.name}); %
    %=======
    %   seasonal
    %=======
    %{
    formats1=['%d',repmat('%f',[1,3]),repmat('%s',[1,5])];
    for ix=1:size(sublist1,1)  %txt�ļ��ĸ��� 
        fid1=fopen([subpath,sublist1(ix).name]);
        a=textscan(fid1,formats1,'headerlines',1,'delimiter','\t','CollectOutput',1);%        
        fclose(fid1);
        dat3=a{1,3};%91*5cell
        aimdat=cell2mat(dat3);%�ܳ�������repmat('%s',[1,5]���ļ����NaN���Զ�ת��Ϊ0
        
    end
    %}
    for ix=1:size(sublist1,1)  %txt�ļ��ĸ��� 
        a=dlmread([subpath,sort_nat_name1{ix}],'\t',1,0);%���������ļ�������
        %a=dlmread([subpath,sublist1(ix).name],'\t',1,0);%���Ͻǵ���ʼ������R=0,C=0 ���õ�91*9��double����  
        %a2=importdata([subpath,sublist1(ix).name],'\t',1);%1��ʾȡ��һ�����£��õ�һ���ṹ����
        %a is tag-lat-lon-ht-���ļ�
        season_mins=nanmean(a(:,5:end));%nanmin(a(:,5:end));
        season_maxs=nanmax(a(:,5:end));
        aout=round([season_mins;season_maxs],0);
        fid1=fopen([matpath,'\',fname1],'a');%��׷��ģʽ��  
        fprintf(fid1,'%s\r\n',sort_nat_name1{ix});%'%s\t'��������ַ������һ��tab�����ø�������ַ����ֿ�
        fprintf(fid1,'%d  %d  %d  %d  %d\r\n',aout');
        fclose(fid1);
        disp([i,ix])  %ȷ���������������ļ���
%         dlmwrite([matpath,'\',fname1],[sublist1(ix).name,'\n'],'-append');%dlmwrite�Զ��ָ��ַ���          
%         dlmwrite([matpath,'\',fname1],aout,'delimiter','\t','-append');        
    end   
    
    %=======
    %   monthly
    %=======    
    for iy=1:size(sublist1,1)  %txt�ļ��ĸ��� 
        b=dlmread([subpath,sort_nat_name2{iy}],'\t',1,0);%���������ļ�������
        %b=dlmread([subpath,sublist2(iy).name],'\t',1,0);%���Ͻǵ���ʼ������R=0,C=0 ���õ�91*9��double����  
        %b is tag-lat-lon-ht-1��~12��
        mon_mins=nanmean(b(:,5:end));%
        mon_maxs=nanmax(b(:,5:end));
        bout=round([mon_mins;mon_maxs],0);
        fid2=fopen([matpath,'\',fname2],'a');%��׷��ģʽ��  
        fprintf(fid2,'%s\r\n',sort_nat_name2{iy});%'%s\t'��������ַ������һ��tab�����ø�������ַ����ֿ�
        fprintf(fid2,[repmat('%d  ',[1,11]),'%d\r\n'],bout');
        fclose(fid2);
        disp([i,iy])  %ȷ���������������ļ���      
    end     
    
end

%%  6.������λ����ȥ��������תд���ݣ��ļ���������
%ֻ��������ļ�������תд���ɣ��������ݵĽ������
clear,clc
%תдΪ��һ��̨վ���ڶ���γ�ȣ������о��ȣ������е���β���������ꡢ�����ﶬ

list=dir('Z:\data\obs\2001-2021pre-t-rh\pro\matfile\2023*');
%�ֶ����list��������name��Ӧ�����ļ��У�folderΪ�ϼ�·��{'Z:\data\obs\2001-2021pre-t-rh\pro\matfile'}
nf=size(list,1); %��ȡ�ļ��еĸ���  
matpath=list(1).folder;  
outpath='Z:\data\obs\2001-2021pre-t-rh\pro\matfile\forGrads\';
%-----------------

%-----------------

for i=1:nf
    %subfoldername=list.name{i,:};%wrong wrong wrong
    subfoldername=list(i).name;
    subpath=strcat(matpath,'\',subfoldername,'\');
    sublist1=dir([subpath,'��λ��from*�ļ�*.txt']); 
    sublist2=dir([subpath,'��λ��from*����*.txt']); 
    %=======
    %   seasonal
    %=======

    for ix=1:size(sublist1,1)  %txt�ļ��ĸ��� 
        a=dlmread([subpath,sublist1(ix).name],'\t',1,0);%���Ͻǵ���ʼ������R=0,C=0 ���õ�91*9��double����  
        %a2=importdata([subpath,sublist1(ix).name],'\t',1);%1��ʾȡ��һ�����£��õ�һ���ṹ����
        %a is tag-lon-lat-ht-���ļ�
        %sublist1(ix).name is 'from0.1to1.5mm��վ���ļ�vis��ֵ2016-2020.txt'
        aout=[a(:,1),a(:,3),a(:,2),a(:,5:end)];%aout is tag-lat-lon-���ļ� 
        aout=aout'; %תдΪ��һ��̨վ���ڶ���γ�ȣ������о��ȣ������е���β���������ꡢ�����ﶬ
        aout(isnan(aout))=-999;%32766; %0406�޸����
        
        %������fmt1���91��̨վ�󲻻���Ŷ
        %fmt1=repmat('%-8d',[1 91]);%Ĭ���Ҷ��룬�Ӹ�-�ű�Ϊ�����
        %fmt2=repmat('%-8.2f',[1 91]);  
        fmt1=[repmat('%-8d',[1 90]),'%-8d\r\n'];%һ��fmt1ĩβ����
        fmt2=[repmat('%-8.2f',[1 90]),'%-8.2f\r\n'];%һ��fmt2ĩβ����   
        
        %grads��������
        xname=sublist1(ix).name;
        headername=xname(4:length(xname)-23);%ȥ���ʿغ�������
        tailname=xname(length(xname)-12:end);
        midname='_ANNseasonsVISthres';
        a_newname=strcat(headername,midname,tailname);
        
        fid1=fopen([outpath,a_newname],'w+');%�Ը���ģʽ��
        %fid1=fopen([outpath,'trans_',sublist1(ix).name],'w+');%�Ը���ģʽ��
        %fid1=fopen([subpath,'trans_',sublist1(ix).name],'w+');%�Ը���ģʽ��,NaN�����滻  
        %fprintf(fid1,[fmt1,fmt2,fmt2,fmt1,fmt1,fmt1,fmt1,'%8d\r\n'],aout');%old
        fprintf(fid1,[fmt1,fmt2,fmt2,fmt1,fmt1,fmt1,fmt1,fmt1],aout');
        fclose(fid1);
        disp([i,ix])  %ȷ���������������ļ���    
    end    
    
    
    %=======
    %   monthly
    %=======    
    for iy=1:size(sublist1,1)  %txt�ļ��ĸ��� 
        b=dlmread([subpath,sublist2(iy).name],'\t',1,0);%���Ͻǵ���ʼ������R=0,C=0 ���õ�91*9��double����  
        %b is tag-lat-lon-ht-1��~12��
        %bout=[b(:,1),b(:,3),b(:,2),b(:,5:end)];%bout is tag-lat-lon-1��12��
        bout=[b(:,1),b(:,3),b(:,2),round(b(:,5:end),0)];%bout is tag-lat-lon-1��12�£�������round(b(:,5:end),0)
        bout=bout'; %תдΪ��һ��̨վ���ڶ���γ�ȣ������о��ȣ������е���β���������ꡢ�����ﶬ
        bout(isnan(bout))=-999;%32766; %0406�޸����
        
       %grads��������
        xname=sublist1(iy).name; %attention sublist1(iy).name
        headername=xname(4:length(xname)-23);%ȥ���ʿغ�������
        tailname=xname(length(xname)-12:end);
        midname='_MonthlyVISthres';
        b_newname=strcat(headername,midname,tailname);
        
        fid2=fopen([outpath,b_newname],'w+');%�Ը���ģʽ��
        %fid2=fopen([outpath,'trans_',sublist2(iy).name],'w+');%�Ը���ģʽ��
        %fid2=fopen([subpath,'trans_',sublist2(iy).name],'w+');%�Ը���ģʽ��,NaN�����滻
        %fprintf(fid2,[fmt1,fmt2,fmt2,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,'%8d\r\n'],bout');%old
        fprintf(fid2,[fmt1,fmt2,fmt2,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1,fmt1],bout');
        fclose(fid2);
        disp([i,iy])  %ȷ���������������ļ���      
    end     
    
end

%%  ANNorMonthly�ֶ�����   Ӣ���ļ���������дctl
%{
%��Ҫ�˽����Ϣ���Լ��鿴�ļ������Ƶ�
pathin_obsgrd='E:\fortran\stn2grd\';%obs����ת�����ƺ��ŵ�·��
pathin_gridgrd='E:\fortran\stn2grd\netgrids\';%�����ת�����ƺ��ŵ�·��
header_gridgrd='grids_';%�����ת�����ƺ���ļ���ͷ
tail_grd='.txt.grd';%ת������ʱֱ����txt�������.grd

%}

clear,clc
%----��дctl���ļ���----
pathin='Z:\data\obs\2001-2021pre-t-rh\pro\matfile\forGrads\';%���е�seasonal.txt��monthly.txt
%�ֶ��޸���������
filelist='monthly.txt'; %monthly.txt  or  seasonal.txt
keywords='monthly';% seasonal   or monthly

fid=fopen([pathin,filelist],'r');
x=textscan(fid,'%s');%x is 1*1 cell,
xx=x{1,1}; %x{1,1}����13*1cell,������13���ļ���
fclose(fid); 

%----����ctl����------
ctlpathin='Z:\STUDY-Storage-box\grads\GDScripts\medianVIS\';
pathout=ctlpathin;

%�鿴aa��ȷ��վ�����ݵ�ctl����Ҫ�޸ĵڶ�������
fid0=fopen([ctlpathin,keywords,'VIS.ctl'],'r');%
a=textscan(fid0,'%s');%x is 1*1 cell,
aa=a{1,1}; %a{1,1}����22*1cell, 
fclose(fid0); 

%�鿴cc��ȷ��������ctl����Ҫ�޸ĵڶ���
fid0=fopen([ctlpathin,keywords,'VISgrid.ctl'],'r');
c=textscan(fid0,'%s');%x is 1*1 cell,
cc=c{1,1}; %c{1,1}����36*1cell,
fclose(fid0); 

% % %========================================
% % %���ݴ�дctl���ļ���xx���ж�վ�����ݺ������ctl���޸�
% % %========================================
%��Ҫ�˽����Ϣ���Լ��鿴�ļ������Ƶȣ�������ݴ��·�����޸�
pathin_obsgrd='Z:\STUDY-Storage-box\grads\GDScripts\medianVIS\grd\';%obs����ת�����ƺ��ŵ�·��
pathin_gridgrd='Z:\STUDY-Storage-box\grads\GDScripts\medianVIS\grd\';%�����ת�����ƺ��ŵ�·��
header_obsgrd='';%obs����ת�����ƺ���ļ�����ͷ
header_gridgrd='grids_';%�����ת�����ƺ���ļ���ͷ
tail_grd='.txt.grd';%ת������ʱֱ����txt�������.grd

for i=1:size(xx,1)
    fname=xx{i};%eg: from0.1to0.4mm_ANNseasonsVISthres2016-2020.txt
    newname=fname(1:length(fname)-4);%ȥβ��from0.1to0.4mm_ANNseasonsVISthres2016-2020
    tmpaa=aa;
    tmpaa{2}=strcat(pathin_obsgrd,header_obsgrd,newname,tail_grd); %����д��ctl�У�ʶ��̨վ���ݵĶ�����
    tmpaa{6}=strcat(pathout,newname,'.map');
    %%%����������ctl���һ�¶����ctl����
    newaa{1}=[tmpaa{1},'  ',tmpaa{2}];
    newaa{2}=[tmpaa{3},'  ',tmpaa{4}];
    newaa{3}=[tmpaa{5},'  ',tmpaa{6}];
    newaa{4}=[tmpaa{7},'  ',tmpaa{8}];
    newaa{5}=[tmpaa{9},'  ',tmpaa{10}];
    newaa{6}=[tmpaa{11},'  ',tmpaa{12},'  ',tmpaa{13},'  ',tmpaa{14},'  ',tmpaa{15}];
    newaa{7}=[tmpaa{16},'  ',tmpaa{17}];
    newaa{8}=[tmpaa{18},'  ',tmpaa{19},'  ',tmpaa{20},'  ',tmpaa{21}];
    newaa{9}=tmpaa{22};    
    %���obs�����ƶ�Ӧ��ctl
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
    %���grid�����ƶ�Ӧ��ctl
    fid3=fopen([pathout,newname,'grid.ctl'],'w+');
    fprintf(fid3,'%s\r\n',newcc{:});
    fclose(fid3); 
    
end