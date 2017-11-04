file_list=dir('.');
n=length(file_list)-2;
for i=1:n
    file=file_list(i+2).name;
    cd(file);
    delete('c*.nii');
    cd ..;
end