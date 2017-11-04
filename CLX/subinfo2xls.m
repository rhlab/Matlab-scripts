file_list=dir('.');
n=length(file_list)-2;
A={'ID' 'Name' 'Age' 'Sex' 'FlipAngle' 'StudyID' 'ImageType' 'TR' 'TE' 'Thickness' 'Slice' 'Rows' 'Columns' 'ScanDate'};
xlswrite('sub_info',A);
for i=1:n
    sub=file_list(i+2).name;
    cd(sub);
    list=dir('.');
    dcm=list(3).name;
    a=dicominfo(dcm) ;
    b=a.PatientName;
    B={sub b.FamilyName a.PatientAge a.PatientSex a.FlipAngle a.StudyID a.MRAcquisitionType a.RepetitionTime a.EchoTime a.SliceThickness a.ImagesInAcquisition a.Rows a.Columns a.StudyDate};
    s=strcat('A',num2str(i+1));
    cd ..
    xlswrite('sub_info',B,'sheet1',s);
end