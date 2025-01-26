function strut_points_act = find_strut_points(struts,PixelIdxList) 

strut_points={};
w = size(struts,1);
l = size(struts,2);
h = size(struts,3);

[a0,b0,c0] = ind2sub([w,l,h],PixelIdxList);
single_conn=zeros(max(a0)+1-(min(a0)-1)+1,max(b0)+1-(min(b0)-1)+1,max(c0)+1-(min(c0)-1)+1);
a1=a0-(min(a0)-1)+1;
b1=b0-(min(b0)-1)+1;
c1=c0-(min(c0)-1)+1;
ind_local=sub2ind([size(single_conn,1),size(single_conn,2),size(single_conn,3)],a1,b1,c1);



single_conn(ind_local)=1;
end_points = (convn(single_conn,ones(3,3,3),'same')==2) & single_conn;
end_points_ind = find(end_points(:));

if ~isempty(end_points_ind)

    w = size(single_conn,1);
    l = size(single_conn,2);
    h = size(single_conn,3);
    [end_points_coord(:,1),end_points_coord(:,2),end_points_coord(:,3)] = ind2sub([w,l,h],end_points_ind);


    nk=1;
    while size(strut_points,2)<size(PixelIdxList,1)

        if nk==1
        kern=single_conn(end_points_coord(1,1)-1:end_points_coord(1,1)+1,end_points_coord(1,2)-1:end_points_coord(1,2)+1,end_points_coord(1,3)-1:end_points_coord(1,3)+1);
        [a,b,c] = ind2sub([size(kern,1),size(kern,2),size(kern,3)],find(kern(:)));
        a=a+end_points_coord(1,1)-1-1;
        b=b+end_points_coord(1,2)-1-1;
        c=c+end_points_coord(1,3)-1-1;
       neigh_local_ind=sub2ind([size(single_conn,1),size(single_conn,2),size(single_conn,3)],a,b,c);
        new_neigh = setdiff(neigh_local_ind,end_points_ind(1));
       strut_points{1}=end_points_ind(1);
       strut_points{2}=new_neigh;
        nk=2;
        else 
          nk=nk+1;
          [new_neigh_coord(1,1),new_neigh_coord(1,2),new_neigh_coord(1,3)]=ind2sub([w,l,h],new_neigh);
         if new_neigh~=end_points_ind(2)
                kern=single_conn(new_neigh_coord(1,1)-1:new_neigh_coord(1,1)+1,new_neigh_coord(1,2)-1:new_neigh_coord(1,2)+1,new_neigh_coord(1,3)-1:new_neigh_coord(1,3)+1);
                [a,b,c] = ind2sub([size(kern,1),size(kern,2),size(kern,3)],find(kern(:)));
                 a=a+new_neigh_coord(1,1)-1-1;
                b=b+new_neigh_coord(1,2)-1-1;
                c=c+new_neigh_coord(1,3)-1-1;
                neigh_local_ind=sub2ind([size(single_conn,1),size(single_conn,2),size(single_conn,3)],a,b,c);               
                new_neigh = setdiff(neigh_local_ind,[strut_points{:}]');
                strut_points{nk}=new_neigh; 
         end
    
        end
    end    
 
    [a,b,c] = ind2sub([w,l,h],[strut_points{:}]');
    [strut_points_act] = sub2ind([size(struts,1),size(struts,2),size(struts,3)],a+(min(a0)-1)-1,b+(min(b0)-1)-1,c+(min(c0)-1)-1);
else
    
    strut_points_act=[];
end
    
end