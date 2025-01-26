function flag = identify_border_struts(strut_points_act,strut_points_thick,struts_jun_bi) 

%% identify all free skeletal branches that touch the boundary of the microstructure domain

    w = size(struts_jun_bi,1);
    l = size(struts_jun_bi,2);
    h = size(struts_jun_bi,3);

    
    [a,b,c] = ind2sub([w,l,h],strut_points_act);
    max_thick = max(strut_points_thick);
    min_mat=cat(2,a-max_thick,b-max_thick,c-max_thick);
    max_mat=cat(2,a+max_thick,b+max_thick,c+max_thick);
    
    if ~isempty(find(min_mat<=1)) || ~isempty(find(max_mat(:,1)>=w)) || ~isempty(find(max_mat(:,2)>=l)) || ~isempty(find(max_mat(:,3)>=h))
        flag =1;
    else
        flag=0;
    end
    
end   
    
    
    