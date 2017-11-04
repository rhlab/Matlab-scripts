file_list=dir('.');
n=length(file_list)-3;
for i=1:n
    a=load_nii(file_list(i+3).name);
    a.img=abs(a.img);
    cd a;
    save_nii(a,file_list(i+3).name);
    cd ..;
    clear a;
end
clear i file_list n;
