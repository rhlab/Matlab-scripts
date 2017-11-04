file_list=dir('.');
n=length(file_list)-2;

for i=1:n
    sub=file_list(i+2).name;
    cd(sub);
    file=dir('.');
    if length(file)>400
        delete('*.txt');
    end
    clear file
    cd ..
end