function [strut_segments,aspect_ratio,max_eccen,area_min,min_thick,len] = find_strut_segments(strut,area,thick,threshold_eccentricity,threshold_aspect_ratio,threshold_cross_area,threshold_change_area)


area_vec=[area{:}]';
thick_vec=thick;


area_min=min(area_vec);
if area_min<threshold_cross_area    %% apply threshold of minimum cross-sectional area
    area_id=find(area_vec<(area_min*threshold_change_area));  %% isolate slender region frmo non-slender (junction) region
    area_id=area_id';
    p=find(diff(area_id)>1);
    area_ind=[area_id(1),area_id(p+1);area_id(p),area_id(end)];
    set_ids={};
    for i=1:size(area_ind,2)
        set_ids{i}=(area_ind(1,i):area_ind(2,i));    
    end
    
    nid=0;
    id=[];
    for i=1:size(set_ids,2)
        
        if ~isempty(intersect(area_vec(set_ids{i}),area_min))
            eq_rad=sqrt(area_vec(set_ids{i})/pi);
            eccen=eq_rad./thick_vec(set_ids{i}); %% calculate eccentricity at each skeletak point
            min_thick=min(thick_vec(set_ids{i}));
            len=length(set_ids{i});
            aspect_ratio=len/(2*min_thick); %% calculate aspect ratio of the strut
            max_eccen=max(eccen);  %% calculate maximum eccentricity
    
            if max_eccen<threshold_eccentricity && aspect_ratio>threshold_aspect_ratio %% apply threshold of 1:eccentricity and  2: aspect ratio 
               nid=nid+1;
              id(nid)=i;
            end
        end
     end
     for i=1:size(set_ids,2)
           size_set(i)=size(set_ids{i},2);   
     end

     [~,maxid]=max(size_set(id));

     if isempty(maxid)
           strut_segments=[];
            
     else
         strut_segments=strut(set_ids{id(maxid)});
     end
     
else
     strut_segments=[];
     aspect_ratio=[];
     max_eccen=[];
     min_thick=[];
     len = [];
end




end