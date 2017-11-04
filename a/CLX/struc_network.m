function [Means stds P_value Dvalue Dvalue_Raw Num net_g1g2 net_all] = struc_network(g1n,g2n,deltas,sparstiyNum,repeatNum,brainRegion)


%==========================================================================
% structural network calculate 
%
% Note：测试版本为默认68个脑区，其他脑区数量请更改相应位置
%
% Means: 随机1000次之后每个稀疏度下每个网络参数（有7个网络参数）的平均值
%
% stds: 随机1000次之后每个稀疏度下每个网络参数（有7个网络参数）的标准差
%
% P_value: 最初的真实值（真正的两个分组）与随机1000次的随机值比，每个稀疏度下每个网络参数的P值
%
% Dvalue: 随机1000次后两组7个网络参数值之差
%
% Dvalue_Raw: 最初正常被试与病人计算7个网络参数值之差
%
% Num: 随机1000次，每次的每个稀疏度下大于真实值之差的个数
%
% net_g1g1: 存放真实值的网络参数
%
% net_all: 存放重复1000次每次的网络参数
%
% reference by 
%
% Structural Insights into Aberrant Topological Patterns of Large-Scale 
%
% Cortical Networks in Alzheimer's Disease
%
%
% Input：g1n number of group1 g1n输入第一组被试人数
%        g2n number of group2 g2n输入第二组被试人数
%        deltas 稀疏度间隙多少
%        sparstiyNum 输入多少个稀疏度 (Smin+deltas*sparstiyNum,例如 25个间隙，sparstiyNum应该输入24)
%        repeatNum 重复次数
%        brainRegion 有多少个脑区
%
% 例子：[Means stds P_value Dvalue Dvalue_Raw Num net_g1g2 net_all] = struc_network(78,90,0.02,24,1000,68)
%
% Writing by Xiaojin Liu 2016/5/2 13:45
% Renew by Meiqi Niu 2016/6/17 15:20
%==========================================================================

%以下开始计算正常分组情况下两组被试脑网络参数
%首先在工作路径下新建一个mat格式文件，命名为all_sub.mat，包括两组被试所有脑区的皮层厚度
load all_sub

MTCG1 = all_sub(1:g1n,:);
MTCG2 = all_sub(g1n+1:g1n+g2n,:);

%计算正常人和病人两个相关矩阵 第一个是fcg1_r，第二个是fcg2_r
Nvar = size(MTCG1,2);
    fcg1_r = zeros(Nvar); fcg1_p = zeros(Nvar);
    for regi = 1:Nvar-1;
        for regj = regi+1:Nvar;
            [R P] = corrcoef(MTCG1(:,regi),MTCG1(:,regj));
            fcg1_r(regi,regj) = R(2); fcg1_p(regi,regj) = P(2);
        end;
    end;
    fcg1_r = fcg1_r + fcg1_r';
    fcg1_p = fcg1_p + fcg1_p';
    
  clear Nvar regi regj R P;
    
 Nvar = size(MTCG2,2);
    fcg2_r = zeros(Nvar); fcg2_p = zeros(Nvar);
    for regi = 1:Nvar-1;
        for regj = regi+1:Nvar;
            [R P] = corrcoef(MTCG2(:,regi),MTCG2(:,regj));
            fcg2_r(regi,regj) = R(2); fcg2_p(regi,regj) = P(2);
        end;
    end;
    fcg2_r = fcg2_r + fcg2_r';
    fcg2_p = fcg2_p + fcg2_p';
    
  clear Nvar regi regj R P;

%以下开始进行网络参数计算的常规步骤
Ws={fcg1_r;fcg2_r};

for i=1:2;
    m=find(Ws{i,1}<0);
    Ws{i,1}(m)=0;
end;
  
  for i=1:2;
    Nm=find(Ws{i,1}<0);
  end;
 
  Sparstiy_Def=cell(2,1);
  Rmatrix=Ws{1,1};
  [Rmax Smin Kmin] = gretna_get_rmax (Rmatrix);
  Sparstiy_Def{1,1}= [Rmax Smin Kmin];
  Rmatrix=Ws{2,1};
  [Rmax Smin Kmin] = gretna_get_rmax (Rmatrix);
  Sparstiy_Def{2,1}= [Rmax Smin Kmin];
  
%s1为最小稀疏度的值，这个由计算得到，deltas为稀疏度的间隙
s1=Smin;
s2=Smin+deltas*sparstiyNum; %多少个稀疏度间隙
n=100;
Thres_type='s';
[net node] = gretna_sw_batch_networkanalysis_weight(Ws, s1, s2, deltas, n, Thres_type);

Dvalue_Raw_Cp=net.Cp(:,1)-net.Cp(:,2);
Dvalue_Raw_Lp=net.Lp(:,1)-net.Lp(:,2);
Dvalue_Raw_gE=net.gE(:,1)-net.gE(:,2);
Dvalue_Raw_locE=net.locE(:,1)-net.locE(:,2);
Dvalue_Raw_Cpratio=net.Cpratio(:,1)-net.Cpratio(:,2);
Dvalue_Raw_Lpratio=net.Lpratio(:,1)-net.Lpratio(:,2);
Dvalue_Raw_Sigma=net.Cpratio(:,1)./net.Lpratio(:,1)-net.Cpratio(:,2)./net.Lpratio(:,2);

for g=1:brainRegion; %注意更改脑区数量，if之后条件语句记得+1！！！！！
    if g<brainRegion+1; % if之后条件语句记得+1！！！！！
Dvalue_Raw_node_gE(:,g)=node.gE(:,1,g)-node.gE(:,2,g);
Dvalue_Raw_node_deg(:,g)=node.deg(:,1,g)-node.deg(:,2,g);
Dvalue_Raw_node_bw(:,g)=node.bw(:,1,g)-node.bw(:,2,g);
    end;
end;



Dvalue_Raw={Dvalue_Raw_Cp,Dvalue_Raw_Lp,Dvalue_Raw_gE,Dvalue_Raw_locE,Dvalue_Raw_Cpratio,Dvalue_Raw_Lpratio,Dvalue_Raw_Sigma,Dvalue_Raw_node_gE,Dvalue_Raw_node_deg,Dvalue_Raw_node_bw}

net_g1g2={net}
node_g1g2={node}

clear net Rmax Smin Kmin Ws;

%下面开始随机分配矩阵，重复1000次计算
for j=1:repeatNum;
    if j<repeatNum+1;
idx=randperm(g1n+g2n);
idx=idx(1:g1n);
MTC1=all_sub(idx,:);
MTC2=all_sub;
MTC2(idx,:)=[];

Nvar = size(MTC1,2);    
    fc1_r = zeros(Nvar); fc1_p = zeros(Nvar);
    for regi = 1:Nvar-1;
        for regj = regi+1:Nvar;                 
            [R P] = corrcoef(MTC1(:,regi),MTC1(:,regj));
            fc1_r(regi,regj) = R(2); fc1_p(regi,regj) = P(2);            
        end;
    end;
    fc1_r = fc1_r + fc1_r';
    fc1_p = fc1_p + fc1_p';
    
  clear Nvar regi regj R P;
    
 Nvar = size(MTC2,2);    
    fc2_r = zeros(Nvar); fc2_p = zeros(Nvar);
    for regi = 1:Nvar-1;
        for regj = regi+1:Nvar;                      
            [R P] = corrcoef(MTC2(:,regi),MTC2(:,regj));
            fc2_r(regi,regj) = R(2); fc2_p(regi,regj) = P(2);            
        end;
    end;
    fc2_r = fc2_r + fc2_r';
    fc2_p = fc2_p + fc2_p';  
    
  clear Nvar regi regj R P;
  
  Ws={fc1_r;fc2_r};
  
  clear fc1_r fc2_r
  
  for i=1:2;
    m=find(Ws{i,1}<0);
    Ws{i,1}(m)=0;
end;
  
  for i=1:2;
    Nm=find(Ws{i,1}<0);
  end;
 
[net node] = gretna_sw_batch_networkanalysis_weight(Ws, s1, s2, deltas, n, Thres_type);

Dvalue_Cp(:,j)=net.Cp(:,1)-net.Cp(:,2);
Dvalue_Lp(:,j)=net.Lp(:,1)-net.Lp(:,2);
Dvalue_gE(:,j)=net.gE(:,1)-net.gE(:,2);
Dvalue_locE(:,j)=net.locE(:,1)-net.locE(:,2);
Dvalue_Cpratio(:,j)=net.Cpratio(:,1)-net.Cpratio(:,2);
Dvalue_Lpratio(:,j)=net.Lpratio(:,1)-net.Lpratio(:,2);
Dvalue_Sigma(:,j)=net.Cpratio(:,1)./net.Lpratio(:,1)-net.Cpratio(:,2)./net.Lpratio(:,2);

for g=1:brainRegion; 
    if g<brainRegion+1; 
Dvalue_node_gE(:,g)=node.gE(:,1,g)-node.gE(:,2,g);
Dvalue_node_deg(:,g)=node.deg(:,1,g)-node.deg(:,2,g);
Dvalue_node_bw(:,g)=node.bw(:,1,g)-node.bw(:,2,g);
    end;
end;

Dvalue_node_gE_all{j,1}=Dvalue_node_gE;
Dvalue_node_deg_all{j,1}=Dvalue_node_deg;
Dvalue_node_bw_all{j,1}=Dvalue_node_bw;

mean_node_gE=mean(Dvalue_node_gE);
mean_node_deg=mean(Dvalue_node_deg);
mean_node_bw=mean(Dvalue_node_bw);

mean_node_gE_all{j,1}=mean_node_gE;
mean_node_deg_all{j,1}=mean_node_deg;
mean_node_bw_all{j,1}=mean_node_bw;

std_node_gE=std(Dvalue_node_gE);
std_node_deg=std(Dvalue_node_deg);
std_node_bw=std(Dvalue_node_bw);

std_node_gE_all{j,1}=Dvalue_node_gE;
std_node_deg_all{j,1}=Dvalue_node_deg;
std_node_bw_all{j,1}=Dvalue_node_bw;

net_all(j,1)={net};
node_all(j,1)={node};


clear net;
clear node;

fprintf('Do not worry, the result is coming coming coming~~~ The %d \n', j)
fprintf('(╯▔︵▔)╯ ≡(▔﹏▔)≡ (￣ε ￣) ╮(╯▽╰)╭ ')




    end;
end;

save net_node;
clear j;

        

for j=1:repeatNum;%生成多少个regions在多少个稀疏度下共多少次的结果(gE)
    if j<repeatNum+1;
    node_a=cell2mat(Dvalue_node_gE_all(j,1));%Dvalue_node_gE_all存放了1000个多个稀疏度下68个脑区的结果，现在分离出来
        for ii=1:brainRegion; 
        if ii<brainRegion+1; 
    str=['region',num2str(ii),'_gE','(:,j)','=node_a(:,ii)'];%region1_gE:多个稀疏度下1000重复的结果，后面以此类推
    eval(str);
    end;
    end;
end;
end;

for j=1:repeatNum;%生成多少个regions在多少个稀疏度下共多少次的结果(deg)
    if j<repeatNum+1;
    node_a=cell2mat(Dvalue_node_deg_all(j,1));
        for ii=1:brainRegion; 
        if ii<brainRegion+1; 
    str=['region',num2str(ii),'_deg','(:,j)','=node_a(:,ii)'];
    eval(str);
    end;
    end;
end;
end;

for j=1:repeatNum;%生成多少个regions在多少个稀疏度下共多少次的结果(bw)
    if j<repeatNum+1;
    node_a=cell2mat(Dvalue_node_bw_all(j,1));
        for ii=1:brainRegion; %注意更改脑区数量
        if ii<brainRegion+1; %if之后条件语句记得+1！！！！！
    str=['region',num2str(ii),'_bw','(:,j)','=node_a(:,ii)'];
    eval(str);
    end;
    end;
end;
end;
        
   
        
        

Dvalue={Dvalue_Cp,Dvalue_Lp,Dvalue_gE,Dvalue_locE,Dvalue_Cpratio,Dvalue_Lpratio,Dvalue_Sigma,Dvalue_node_gE,Dvalue_node_deg,Dvalue_node_bw};

mean_Cp=mean(Dvalue_Cp');
mean_Lp=mean(Dvalue_Lp');
mean_gE=mean(Dvalue_gE');
mean_locE=mean(Dvalue_locE');
mean_Cpratio=mean(Dvalue_Cpratio');
mean_Lpratio=mean(Dvalue_Lpratio');
mean_Sigma=mean(Dvalue_Sigma');




Means={mean_Cp,mean_Lp,mean_gE,mean_locE,mean_Cpratio,mean_Lpratio,mean_Sigma};

std_Cp=std(Dvalue_Cp');
std_Lp=std(Dvalue_Lp');
std_gE=std(Dvalue_gE');
std_locE=std(Dvalue_locE');
std_Cpratio=std(Dvalue_Cpratio');
std_Lpratio=std(Dvalue_Lpratio');
std_Sigma=std(Dvalue_Sigma');



stds={std_Cp,std_Lp,std_gE,std_locE,std_Cpratio,std_Lpratio,std_Sigma};

for k=1:sparstiyNum+1;
    if k<sparstiyNum+2;
Num_Cp=find(Dvalue_Raw_Cp(k,1)<Dvalue_Cp(k,:));
Num_Cp_total(k,1)=length(Num_Cp);
Num_Lp=find(Dvalue_Raw_Lp(k,1)<Dvalue_Lp(k,:));
Num_Lp_total(k,1)=length(Num_Lp);
Num_gE=find(Dvalue_Raw_gE(k,1)<Dvalue_gE(k,:));
Num_gE_total(k,1)=length(Num_gE);
Num_locE=find(Dvalue_Raw_locE(k,1)<Dvalue_locE(k,:));
Num_locE_total(k,1)=length(Num_locE);
Num_Cpratio=find(Dvalue_Raw_Cpratio(k,1)<Dvalue_Cpratio(k,:));
Num_Cpratio_total(k,1)=length(Num_Cpratio);
Num_Lpratio=find(Dvalue_Raw_Lpratio(k,1)<Dvalue_Lpratio(k,:));
Num_Lpratio_total(k,1)=length(Num_Lpratio);
Num_Sigma=find(Dvalue_Raw_Sigma(k,1)<Dvalue_Sigma(k,:));
Num_Sigma_total(k,1)=length(Num_Sigma);
    end;
end;

       
for ii=1:brainRegion; 
        if ii<brainRegion+1; 
    str=['node_gE=','region',num2str(ii),'_gE'];%node_gE暂时存放region**_gE矩阵
    eval(str);
    for k=1:sparstiyNum+1;
    if k<sparstiyNum+2;
    Num_node_gE=find(Dvalue_Raw_node_gE(k,ii)<node_gE(k,:));
    Num_node_gE_total(k,1)=length(Num_node_gE);
    str=['Num_node_gE_total_region',num2str(ii),'=Num_node_gE_total'];
    eval(str);
    end;
    end;
        end;
end;

for ii=1:brainRegion;
        if ii<brainRegion+1; 
    str=['node_deg=','region',num2str(ii),'_deg'];%node_gE暂时存放region**_gE矩阵
    eval(str);
    for k=1:sparstiyNum+1;
    if k<sparstiyNum+2;
    Num_node_deg=find(Dvalue_Raw_node_deg(k,ii)<node_deg(k,:));
    Num_node_deg_total(k,1)=length(Num_node_deg);
    str=['Num_node_deg_total_region',num2str(ii),'=Num_node_deg_total'];
    eval(str);
    end;
    end;
        end;
end;

for ii=1:brainRegion;
        if ii<brainRegion+1; 
    str=['node_bw=','region',num2str(ii),'_bw'];%node_gE暂时存放region**_gE矩阵，找到小于标准的两组差的个数
    eval(str);
    for k=1:sparstiyNum+1;
    if k<sparstiyNum+2;
    Num_node_bw=find(Dvalue_Raw_node_bw(k,ii)<node_bw(k,:));
    Num_node_bw_total(k,1)=length(Num_node_bw);
    str=['Num_node_bw_total_region',num2str(ii),'=Num_node_bw_total'];
    eval(str);
    end;
    end;
        end;
end;
    

Num={Num_Cp_total,Num_Lp_total,Num_gE_total,Num_locE_total,Num_Cpratio_total,Num_Lpratio_total,Num_Sigma_total}

P_Cp=Num_Cp_total./repeatNum;
P_Lp=Num_Lp_total./repeatNum;
P_gE=Num_gE_total./repeatNum;
P_locE=Num_locE_total./repeatNum;
P_Cpratio=Num_Cpratio_total./repeatNum;
P_Lpratio=Num_Lpratio_total./repeatNum;
P_Sigma=Num_Sigma_total./repeatNum;

for ii=1:brainRegion;
        if ii<brainRegion+1; 
    str=['PgE_region',num2str(ii),'=Num_node_gE_total_region',num2str(ii),'./1000'];%PgE_region1:每个稀疏度下每个脑区的P值
    eval(str);
    end;
end;

for ii=1:brainRegion; 
        if ii<brainRegion+1; 
    str=['P_deg_region',num2str(ii),'=Num_node_deg_total_region',num2str(ii),'./1000'];
    eval(str);
    end;
    end;
    
for ii=1:brainRegion;
        if ii<brainRegion+1; 
    str=['Pbw_region',num2str(ii),'=Num_node_bw_total_region',num2str(ii),'./1000'];
    eval(str);
    end;
end;  

P_Cp=P_Cp';
P_Lp=P_Lp';
P_gE=P_gE';
P_locE=P_locE';
P_Cpratio=P_Cpratio';
P_Lpratio=P_Lpratio';
P_Sigma=P_Sigma';



P_value={P_Cp,P_Lp,P_gE,P_locE,P_Cpratio,P_Lpratio,P_Sigma}

save structural_network



   