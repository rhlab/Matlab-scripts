% calculate the centroid of ROI in MNI space
% Input: 
%     fn - nii file name
% Output:
%     centroidMNI: centroid in MNI space
%     c: centroid in voxel space
%     origin: the image origin
% Required: 
%     NIFTI toolbox by Jimmy Shen
%% 20151026

function [centroidMNI,c,origin]=ROIcentroid(fn)
  nii=load_nii(fn);
  pixdim=nii.hdr.dime.pixdim(2:4);
  [l,m,n]=ind2sub(size(nii.img),find(nii.img~=0));
  
  c=mean([l m n]);
  
  origin=[nii.hdr.hist.srow_x(4) nii.hdr.hist.srow_y(4) nii.hdr.hist.srow_z(4)];
  centroidMNI=round(c.*pixdim+origin - pixdim);
  
  c=round(c);
end

%% OLD code
% % 
% % clc
% % clear
% % cd H:\tp
% % !dir /b *.nii > list.txt
% % file_list=textread('list.txt','%s');
% % for f=1:length(file_list);
% % files=file_list{f};
% % vol=spm_vol(files);
% % img=spm_read_vols(vol);
% % origin=vol.mat(1:3,4)';
% % [m n p] = size(img);
% % coordinates = zeros(1,3);
% % z = 1;
% % for i = 1:m
% %     for j = 1:n
% %         for k = 1:p
% %             if img(i,j,k)~= 0                
% %                coordinates(z,1) = i;
% %                coordinates(z,2) = j;
% %                coordinates(z,3) = k;
% %                z = z + 1;
% %           end
% %         end
% %     end
% % end
% % center_coord=mean(coordinates);
% % center_coord_MNI(f,:)=center_coord+origin;
% % end
% % save center_coord_MNI center_coord_MNI file_list
% % 
