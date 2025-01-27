% Segment slender regions in 3D micorstructures
%==========================================================================
% AUTHOR        Vinit Vijay Deshpande
% CONTACT       vinit-vijay.deshpande@h-da.de
% INSTITUTION   University of Applied Sciences, Darmstadt

% INPUTS:

%       IMAGE:   3D binary image of the microstructure

%       Threshold values
% 
%           threshold_length: threshold value for pruning based on length
%
%           threshold_loop_ratio : threshold value for pruning based on loop ratio
%
%           threshold_eccentricity: threshold value for pruning based on eccentricity
%           of the cross-sectional area 
%
%           threshold_aspect_ratio: threshold value for pruning based on aspect ratio
%
%           threshold_cross_area: threshold value for pruning based on
%           minimum crosss-sectional area
%
%           threshold_change_area: threshold value for calculating the truncated
%           skeletal branch that separates the slender region fom the junction region
% 

% OUTPUTS:

%       strut_segments: a cell array of all the remaining skeletal branches after
%       pruning
%       strut_pixels: a cell array of all the voxels of the struts segmented by the
%       algorithm

%%%%%%%%%%%%%%%%%%% Select material microstructure %%%%%%%%%%%%%%%%%%%%%%%%
prompt = "Select microstructure (1 or 2): 1. FOAM , 2.IWP";
choice = input(prompt);
switch choice
    case 1
        load FOAM_MICROSTRUCTURE.mat  %% load binary image of the foam micorstructure
    case 2
        load IWP_MICROSTRUCTURE.mat  %% load binary image of the IWP micorstructure
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%Input threshold values%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threshold_length = 3;
threshold_loop_ratio = 0.2;
threshold_eccentricity = 10;
threshold_aspect_ratio = 0.1;
if choice==1
    threshold_cross_area = 4000;
    threshold_change_area = 5;
else if choice==2
    threshold_cross_area = 550;
    threshold_change_area = 1.3;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ori=IMAGE;
ori_s=smooth3(ori,"gaussian",5); %% smoothen the image 
ori=logical(ori_s);

struts = bwskel(ori);    %% create 3D image of skeleton of the  microstructure
struts_dist = bwdist(~ori); %% distance transform
thick = struts.*struts_dist;  %% create a 3D matrix showing thickness values at skeleton points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%to visualize the skeleton and the original image%%%%%%%%%%%%%%

%% select a small portion of the microstructure and the skeleton

% if choice ==1
% %%%%%% FOAM MICROSTRUCTURE
% ori_l = ori(200:400,200:400,200:400);  %% a small region of the microstructure
% struts_l = struts(200:400,200:400,200:400);
% else if choice ==2
% %%%%%% IWP MICROSTRUCTURE
% ori_l = ori(1:100,1:100,1:100);  %% a small region of the microstructure
% struts_l = struts(1:100,1:100,1:100);
% end
% end
% % 
% figure(1);
% col=[.7 .7 .8];
% hiso = patch(isosurface(ori_l,0),'FaceColor',col,'EdgeColor','none');
% axis equal;axis off;
% lighting phong;
% isonormals(ori_l,hiso);
% alpha(0.7);
% set(gca,'DataAspectRatio',[1 1 1])
% camlight;
% hold on;
% w=size(struts_l,1);
% l=size(struts_l,2);
% h=size(struts_l,3);
% [x,y,z]=ind2sub([w,l,h],find(struts_l(:)));
% plot3(y,x,z,'square','Markersize',1,'MarkerFaceColor','k','Color','k');            
% set(gcf,'Color','white');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%determine junction points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

struts_jun =  convn(struts,ones(3,3,3),'same').*struts;
struts_jun_bi=zeros(size(struts_jun));
struts_jun_bi(struts_jun(:)>3)=1;  %% 3D image with all junction points
struts_nojun = ~struts_jun_bi.*struts; %% 3D skeleton image without junction points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%create a cell array of skeletal branches and their thickness%%%%

struts_CC = bwconncomp(struts_nojun);  %%  determine connected components in the skeleton image
strut_points_act={};  %% cell array of skeletal branches
strut_points_thick={};  %% cell array of thickness of skeletal branches
for i=1:size(struts_CC.PixelIdxList,2)
    if size(struts_CC.PixelIdxList{i},1)==1
        strut_points_act{i}=struts_CC.PixelIdxList{i};
    else
        strut_points_act{i} = find_strut_points(struts,struts_CC.PixelIdxList{i});   %% arrange the skeletal points in order
    end
        strut_points_thick{i} = thick(strut_points_act{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% identify and prune free branches %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

strut_points_act_keep={}; %% cell array of remaining skeletal branches after pruning
strut_points_thick_keep={}; %% cell array of thickness of the remaining branches
w = size(struts_jun_bi,1);
l = size(struts_jun_bi,2);
h = size(struts_jun_bi,3);
struts_jun_bi_pad=zeros(w+2,l+2,h+2);
struts_jun_bi_pad(2:end-1,2:end-1,2:end-1)=struts_jun_bi;
n0=0;
n0b=0;
n0c=0;
for i=1:size(strut_points_act,2)
    % i
    flag1= identify_corner_struts(w,l,h,strut_points_act{i},struts_jun_bi_pad);  %% identify free skeletal branches
    
    if flag1==1
        n0c=n0c+1;
        flag2 = identify_border_struts(strut_points_act{i},strut_points_thick{i},struts_jun_bi);   %% identify free skeletal branches at the boundary
        
        if flag2==1
            
            n0=n0+1;
            n0b=n0b+1;
            strut_points_act_keep{n0}=strut_points_act{i};   %% store free skeletal branches that are at the boundary
            strut_points_thick_keep{n0}=strut_points_thick{i};
        end
    else
         n0=n0+1;
         strut_points_act_keep{n0}=strut_points_act{i};   %% store skeletal branches that are not free (refer to manuscript for the definition of a free branch)
         strut_points_thick_keep{n0}=strut_points_thick{i};
   
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%calculate length of the skeletal branches %%%%%%%%%%%%%%%%


size_struts=zeros(size(strut_points_act_keep,2),1);
for i=1:size(strut_points_act_keep,2)
    
    if ~isempty(strut_points_act_keep{i})
        
        size_struts(i)=size(strut_points_act_keep{i},1);
    else 
        size_struts(i)=0;
    end
           
end

% figure(1)
% his=histogram(size_struts,'BinWidth',1);  %% histogram plot of branch length

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%pruning based on length of the skeletal branch%%%%%%%%%%%%%
strut_points_act_long={};
strut_points_thick_long={};
n1=0;
for i=1:size(strut_points_act_keep,2)
    
    if size(strut_points_act_keep{i},1)>threshold_length
      
        n1=n1+1;
        strut_points_act_long{n1}=strut_points_act_keep{i};
        strut_points_thick_long{n1}=strut_points_thick_keep{i};
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%pruning based on loop ratio of the skeletal branch%%%%%%%%%%%%


loop_ratio=zeros(size(strut_points_act_long,2),1);
n2=0;
strut_points_act_straight={};
strut_points_thick_straight={};
for i=1:size(strut_points_act_long,2)
    
    [d,e,f]=ind2sub(size(ori),strut_points_act_long{i});
    [end_point_dist,arc_length]=loop_detection([e,d,f]);
    loop_ratio(i)=end_point_dist/arc_length;
    
    if loop_ratio(i)>threshold_loop_ratio
       n2=n2+1;
       strut_points_act_straight{n2}=strut_points_act_long{i};
       strut_points_thick_straight{n2}=strut_points_thick_long{i};
    end          
end 

% figure(2)
% his1=histogram(loop_ratio,'BinWidth',0.02);  %% histogram plot of loop ratio

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%calculate cross-sectional area (cross_area)of the skeletal branch%%%%%%%

% % this section is resource intensive for FOAM_MICROSTRUCTURE. Please load
% % cross_area.mat and strut_segments.mat if you dont want to run this
% % section (only for FOAM_MICROSTRUCTURE)

cross_area={}; %% cell array of cross-sectional area of the  skeletal branches 
strut_segments={};  %% cell array of truncated skeletal branches
ori_pad=padarray(ori,[200 200 200],0,'both');
size_ori=size(ori);

parfor i=1:size(strut_points_act_straight,2)
    [x,y,z]=ind2sub(size_ori,strut_points_act_straight{i});    
    x=x+200;
    y=y+200;
    z=z+200;
    ori_local=ori_pad(min(x)-200:max(x)+200,min(y)-200:max(y)+200,min(z)-200:max(z)+200);
    x=x-(min(x)-200)+1;
    y=y-(min(y)-200)+1;
    z=z-(min(z)-200)+1;
    [id_x,id_y,id_z]=ind2sub(size(ori_local),1:length(ori_local(:)));
    
    for j=1:length(strut_points_act_straight{i})

        normal=find_normal(strut_points_act_straight{i},x,y,z,j); %% calculate tangent vector to the strut segment at that skeletal point 
        resi=normal(1)*(id_x-x(j))+normal(2)*(id_y-y(j))+normal(3)*(id_z-z(j));
        cross_x=id_x(resi<1 & resi>-1);
        cross_y=id_y(resi<1 & resi>-1);
        cross_z=id_z(resi<1 & resi>-1);
        cross_ind=sub2ind(size(ori_local),cross_x,cross_y, cross_z);
        pixel_value=ori_local(cross_ind);
        white_pixels=find(pixel_value==1);
        trail_image=zeros(size(ori_local));
        trail_image(cross_ind(white_pixels))=1;
         
        cross_x_white=cross_x(white_pixels);
        cross_y_white=cross_y(white_pixels);
        cross_z_white=cross_z(white_pixels);
        
        origin_ind=sub2ind(size(ori_local),x(j),y(j),z(j));
        origin_row=find(cross_ind==origin_ind);
        white_pixles_row=find(white_pixels==origin_row);
        
        [coeff, score]=pca([cross_x_white',cross_y_white',cross_z_white']);
        new_coord_trans=score-min(score)+1;
        max_new_coord_trans=round(max(new_coord_trans));
        new_cross_image=zeros(max_new_coord_trans(1)+10,max_new_coord_trans(2)+10);
        new_coords_ind=sub2ind(size(new_cross_image),round(new_coord_trans(:,1)),round(new_coord_trans(:,2)));
        new_cross_image(new_coords_ind)=1;
        cross_image_filled = (conv2(new_cross_image,ones(3,3),'same')>=1);
        [L, blobCount] = bwlabel(cross_image_filled);
        stats=regionprops(L,'Centroid','Circularity','Eccentricity');
        st_pt_1=round(new_coord_trans(white_pixles_row,1));
        st_pt_2=round(new_coord_trans(white_pixles_row,2));
        st_point=sub2ind(size(cross_image_filled),st_pt_1,st_pt_2);
        blobCount_st=L(st_point);
        if blobCount_st~=0
            cross_area{i}{j}=sum(L(:)==blobCount_st);
        else 
            cross_area{i}{j}=[];
        end       
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%pruning based on aspect ratio, eccentricity and %%%%%%%%%%%%%%%%
%%%%%%%%%%cross-sectional area of strut of the skeletal branch %%%%%%%%%%%%


strut_segments={};  %% a cell array of all remaining skeletal branches after applying all pruning stratgies
strut_aspect_ratio={}; %% a cell array of aspect ratio of struts correspoding to all remaining sskeletal branches
nsg=0;
max_eccen={}; %% a cell array of maximum eccentricity of struts corresponding to all remaining sskeletal branches
area_min={};  %% a cell array of minimum cross-sectional area of struts corresponding to all remaining skeletal branches
min_thick={};  %% a cell array of minimum thickness of struts corresponding to all remaining skeletal branches
len_strut_segments={}; %% a cell array of length of struts corresponding to all remaining skeletal branches

parfor i=1:size(strut_points_act_straight,2)  
    
    [strut_segments{i},strut_aspect_ratio{i},max_eccen{i},area_min{i},min_thick{i},len_strut_segments{i}] = find_strut_segments(strut_points_act_straight{i},cross_area{i},strut_points_thick_straight{i},threshold_eccentricity,threshold_aspect_ratio,threshold_cross_area,threshold_change_area);  

    if ~isempty(strut_segments{i})
        nsg=nsg+1;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% plot histogram of minimum cross-sectional area of the struts%%%
size_area_min=[];
ns=0;
for i=1:size(area_min,2)
    if ~isempty(strut_segments{i})
        ns=ns+1;
        size_area_min(ns)=area_min{i};
    end
end
figure(3)
his2=histogram(size_area_min,'BinWidth',100);  %% histogram plot of min cross area

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% plot histogram of maximum eccentricity ofcross-sectional area %%% 
% %%%%%%%%%%%%%%%%%%of the struts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

size_max_eccen=[];
ns=0;
for i=1:size(max_eccen,2)
    if ~isempty(strut_segments{i})
        ns=ns+1;
        size_max_eccen(ns)=max_eccen{i};
    end
end
figure(4)
his3=histogram(size_max_eccen,'BinWidth',0.5); %% histogram plot of max eccentricity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% plot histogram of aspect ratio of the struts %%%%%%%%%%%%%%%%%%%% 

size_strut_aspect_ratio=[];
ns=0;
for i=1:size(strut_aspect_ratio,2)
    if ~isempty(strut_segments{i})
        ns=ns+1;
        size_strut_aspect_ratio(ns)=strut_aspect_ratio{i};
    end
end
figure(5)
his6=histogram(size_strut_aspect_ratio,'Normalization','Probability','BinWidth',0.1); %% histogram plot of aspect ratio
ylim([0 0.35])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%% to visualize a particluar skeletal branch in micorstructure %%%


[xcj,ycj,zcj]=ind2sub([size(ori,1),size(ori,2),size(ori,3)],strut_segments{6}); %% strut segment no. 6

extend=35;
ori_pad=padarray(ori,[extend extend extend],0,'both');
xcj=xcj+extend;
ycj=ycj+extend;
zcj=zcj+extend;

tri_ori=ori_pad(min(xcj)-extend:max(xcj)+extend,min(ycj)-extend:max(ycj)+extend,min(zcj)-extend:max(zcj)+extend);
xcj=xcj-(min(xcj)-extend)+1;
ycj=ycj-(min(ycj)-extend)+1;
zcj=zcj-(min(zcj)-extend)+1;

figure(6)
col=[.7 .7 .8];
hiso = patch(isosurface(tri_ori,0),'FaceColor',col,'EdgeColor','none');
axis equal;axis off;
lighting phong;
isonormals(tri_ori,hiso);
alpha(0.4);
set(gca,'DataAspectRatio',[1 1 1])
camlight;
hold on
plot3(ycj,xcj,zcj,'square','Markersize',2,'MarkerFaceColor','k','Color','k');    
set(gcf,'Color','white'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% calculate strut domain (volume) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % % this section is resource intensive for FOAM_MICROSTRUCTURE. Please
% %load strut_pixels.mat if you dont want to run this section (only for FOAM_MICROSTRUCTURE)

strut_pixels={};  %% a cell array to store all voxels that belong to each strut 
ori_pad=padarray(ori,[200 200 200],0,'both');
size_ori=size(ori);

for i=1:size(strut_segments,2)

    if size(strut_segments{i},1)>2

        [x,y,z]=ind2sub(size_ori,strut_segments{i});

        x=x+200;
        y=y+200;
        z=z+200;
        minx=min(x);
        miny=min(y);
        minz=min(z);
        maxx=max(x);
        maxy=max(y);
        maxz=max(z);
        ori_local=ori_pad(minx-200:maxx+200,miny-200:maxy+200,minz-200:maxz+200);   
        x=x-(minx-200)+1;
        y=y-(miny-200)+1;
        z=z-(minz-200)+1;
        [id_x,id_y,id_z]=ind2sub(size(ori_local),1:length(ori_local(:)));
    
        normal1=find_normal(strut_segments{i},x,y,z,1);     %% calculate tangent vector to the strut segment at its first skeletal point    
        resi1=normal1(1)*(id_x-x(1))+normal1(2)*(id_y-y(1))+normal1(3)*(id_z-z(1));
        cross_x1=id_x(resi1<1 & resi1>-1);
        cross_y1=id_y(resi1<1 & resi1>-1);
        cross_z1=id_z(resi1<1 & resi1>-1);
        cross_ind1=sub2ind(size(ori_local),cross_x1,cross_y1, cross_z1);
        ori_local(cross_ind1)=0;
        
        normalend=find_normal(strut_segments{i},x,y,z,length(strut_segments{i}));     %% calculate tangent vector to the strut segment at its last skeletal point   
        resiend=normalend(1)*(id_x-x(end))+normalend(2)*(id_y-y(end))+normalend(3)*(id_z-z(end));
        cross_xend=id_x(resiend<6 & resiend>0);
        cross_yend=id_y(resiend<6 & resiend>0);
        cross_zend=id_z(resiend<6 & resiend>0);
        cross_indend=sub2ind(size(ori_local),cross_xend,cross_yend, cross_zend);
        ori_local(cross_indend)=0;
        
        
        [L, blobCount] = bwlabeln(ori_local);
        jn=round(length(strut_segments{i})/2);
        st_point=sub2ind(size(ori_local),x(jn),y(jn),z(jn));
        blobCount_st=L(st_point); 
        ori_local_onlystrut=zeros(size(ori_local));
        ori_local_onlystrut(L(:)==blobCount_st)=1;
        ori_struts_whole=zeros(size(ori_pad));

        for j=1:length(strut_segments{i})

            [strut_seg_coord1,strut_seg_coord2,strut_seg_coord3]=ind2sub(size(ori),strut_segments{i}(j));
            strut_seg_coord1=strut_seg_coord1+200;
            strut_seg_coord2=strut_seg_coord2+200;
            strut_seg_coord3=strut_seg_coord3+200;
            xminus=strut_seg_coord1-round(ff*(strut_points_thick_straight{i}(j)+2));
            xplus=strut_seg_coord1+round(ff*(strut_points_thick_straight{i}(j)+2));
            yminus=strut_seg_coord2-round(ff*(strut_points_thick_straight{i}(j)+2));
            yplus=strut_seg_coord2+round(ff*(strut_points_thick_straight{i}(j)+2));
            zminus=strut_seg_coord3-round(ff*(strut_points_thick_straight{i}(j)+2));
            zplus=strut_seg_coord3+round(ff*(strut_points_thick_straight{i}(j)+2));
            [columnsInImage rowsInImage pagesInImage] = meshgrid (1:(xplus-xminus+1),1:(yplus-yminus+1),1:(zplus-zminus+1));
            xc=strut_seg_coord1-xminus+1;
            yc=strut_seg_coord2-yminus+1;
            zc=strut_seg_coord3-zminus+1;
            sphereVoxels = (rowsInImage - yc).^2 + (columnsInImage - xc).^2 + (pagesInImage - zc).^2 <=    (ff*(strut_points_thick_straight{i}(j)+2)).^2;

            ori_struts_whole(xminus:xplus,yminus:yplus,zminus:zplus)=1*(ori_pad(xminus:xplus,yminus:yplus,zminus:zplus) + sphereVoxels);

        end
        ori_struts_whole(ori_struts_whole(:)>0)=1;
        ori_struts_whole_local=ori_struts_whole(minx-200:maxx+200,miny-200:maxy+200,minz-200:maxz+200);

        final = ori_struts_whole_local & ori_local_onlystrut;
        finald=double(final);
        strut_ind=find(final(:)==1);
        [x_local,y_local,z_local]=ind2sub(size(ori_local),strut_ind);
        x_global=x_local+(minx-200)-1-200;
        y_global=y_local+(miny-200)-1-200;
        z_global=z_local+(minz-200)-1-200;
        strut_pixels{i}=cat(2,x_global,y_global,z_global);
             
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% to plot a particular strut domain in the microstructuree %%%%%%%

if choice == 1
id=6;  %% ID of the chosen strut in cell array : strut_segments

[xcj,ycj,zcj]=ind2sub([size(ori,1),size(ori,2),size(ori,3)],strut_points_act_straight{id});
extend=35;
ori_pad=padarray(ori,[extend extend extend],0,'both');
xcj=xcj+extend;
ycj=ycj+extend;
zcj=zcj+extend;
mxcj=(min(xcj));
mycj=(min(ycj));
mzcj=(min(zcj));
maxcj=(max(xcj));
maycj=(max(ycj));
mazcj=(max(zcj));
tri_ori=ori_pad(mxcj-extend:maxcj+extend,mycj-extend:maycj+extend,mzcj-extend:mazcj+extend);
xcj=xcj-(mxcj-extend)+1;
ycj=ycj-(mycj-extend)+1;
zcj=zcj-(mzcj-extend)+1;


figure(7)
col=[.7 .7 .8];
hiso = patch(isosurface(tri_ori,0),'FaceColor',col,'EdgeColor','none');
axis equal;axis off;
xlabel('X axis')
ylabel('Y axis')
zlabel('Z axis')
lighting phong;
isonormals(tri_ori,hiso);
alpha(0.1);
set(gca,'DataAspectRatio',[1 1 1])
camlight;
hold on


strut_region = zeros(size(ori));
strut_region_ind=sub2ind([size(ori,1),size(ori,2),size(ori,3)],strut_pixels{id}(:,1),strut_pixels{id}(:,2),strut_pixels{id}(:,3));
strut_region(strut_region_ind)=1;
strut_region_pad=padarray(strut_region,[extend extend extend],0,'both');
strut_region_local=strut_region_pad(mxcj-extend:maxcj+extend,mycj-extend:maycj+extend,mzcj-extend:mazcj+extend);

figure(7)
hold on
hiso_sr = patch(isosurface(strut_region_local,0),'FaceColor',[.6 .6 .6],'EdgeColor','none');
isonormals(strut_region_local,hiso_sr);
alpha(0.3);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if choice ==2
% 
% ids=[];
% nid=0;
% strut_region = zeros(size(ori));
% for i=1:size(strut_segments,2)
% 
%     if ~isempty(strut_segments{i})
%         nid=nid+1;
%         ids(nid)=i;       
%         strut_region_ind=sub2ind([size(ori,1),size(ori,2),size(ori,3)],strut_pixels{i}(:,1),strut_pixels{i}(:,2),strut_pixels{i}(:,3));
%         strut_region(strut_region_ind)=1;
%     end
% end
% 
% figure(7)
% col=[.7 .7 .8];
% hiso = patch(isosurface(ori,0),'FaceColor',col,'EdgeColor','none');
% axis equal;axis off;
% xlabel('X axis')
% ylabel('Y axis')
% zlabel('Z axis')
% lighting phong;
% isonormals(ori,hiso);
% alpha(0.1);
% set(gca,'DataAspectRatio',[1 1 1])
% camlight;
% hold on
% 
% figure(7)
% hold on
% hiso_sr = patch(isosurface(strut_region,0),'FaceColor',[.6 .6 .6],'EdgeColor','none');
% isonormals(strut_region,hiso_sr);
% alpha(0.3);
% 
% end
%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

