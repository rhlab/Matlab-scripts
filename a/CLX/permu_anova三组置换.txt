function [pvalue,sig]=permu_anova1(a,b,c,sim,alpha)
%原始输入向量a,b,c,分别代表三个样本，不管是否等组。group为分组向量,sim为置换次数，alpha为显著水平
%该置换检验使用的统计参数为one-way anova的F值
%注意，缺陷在于randperm指令随机生成sim个数列（一般为1000），可能不能保证没有重复
%但permutation检验可能允许该重复，因为其组合数目远大于1000，重复的可能性较低
%pvalue: calculated p value. sig: sig=1, there is a signficant
%difference between a and b; sig=0, no difference. Sim:times of
%randomizations. Alpha:significance level, e.g.,0.05

d=[a;b;c];s1=length(a);s2=length(b);
n=length(d);m=ones(n,1);index1=m(1:s1);index2=2*m(s1+1:s1+s2);
index3=3*m(s1+s2+1:n);group=[index1;index2;index3];
[p anovatab]=anova1(d,group,'off'); 
F0=anovatab{2,5};  %求得原始F0值
clear p anovatab

F=zeros(1,sim);
for ii=1:sim
    ind=randperm(n);                    %随机产生从1到n的排列数值，为一行向量
    anew=d(ind(1:s1));                  %在总样本中取上面生成的随机数列的前s1个（为样本a的元素数目），定义为随机新样本a
    bnew=d(ind(s1+1:s1+s2));
    cnew=d(ind(s1+s2+1:n));
    dnew=[anew;bnew;cnew];
    [p anovatab]=anova1(dnew,group,'off');   %计算随机产生的三个新样本的anova1检验
    F(ii)=anovatab{2,5};                  %求得新的F值
    clear anew bnew ind p anovatab
end

%计算大于等于原始F0值的随机生成的F的数量 
    greaternumber=length(find(F>F0)); 
    pvalue=(greaternumber)/(sim);               %计算显著性水平，用上面数量除以置换次数

if (pvalue<=alpha) sig=1;
else sig=0;
end