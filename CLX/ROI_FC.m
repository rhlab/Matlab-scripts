%individual matrix to cell
for i=1:53
    a=sprintf('ROICorrelation_FisherZ_b%03i',i)
    load(a)
    b{i,1}=ROICorrelation_FisherZ
    clear ROICorrelation_FisherZ
end

%cell to matrix
function [fc_mat]=3fc_ext(fc,s,r)
% fc: fc matrix
% s: number of subjects
% r: number of ROIszaoqi z

for n=1:s   %subjects loop
    a=1;
    for i=1:r-1   %row loop
        for j=1:r-i    %column loop
            fc_mat(n,a)=fc{n,1}(i,i+j);
            a=a+1;
        end
    end
end

%dpabi Ancoval
[F P]=y_ancova1(DependentVariable,GroupLabel,Covariates)

F_FC=ones(21,1)
p_FC=ones(21,1)
for i=1:21
    [F P]=y_ancova1(FC(:,i),label,covs)
    F_FC(i,1)=F
    p_FC(i,1)=P
    clear F P
end
    
    
