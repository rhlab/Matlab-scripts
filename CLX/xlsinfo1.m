function xlsinfo(n)
A={'ID' 'Name' 'Age' 'Sex' 'FlipAngle' 'StudyID' 'ImageType' 'TR' 'TE' 'Thickness' 'Slice' 'Rows' 'Columns' 'ScanDate'};
xlswrite('sub_data',A);
for i=1:n
    sub=sprintf('sub%02i',i);
    cd(sub)
    list=dir('.')
    dcm=list(3).name
    a=dicominfo(dcm) ;
    b=a.PatientName;
    B={sub b.FamilyName a.PatientAge a.PatientSex a.FlipAngle a.StudyID a.MRAcquisitionType a.RepetitionTime a.EchoTime a.SliceThickness a.ImagesInAcquisition a.Rows a.Columns a.StudyDate};
    s=strcat('A',num2str(i+1));
    cd ..
    xlswrite('sub_data',B,'sheet1',s);
end
