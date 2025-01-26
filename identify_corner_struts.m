function flag = identify_corner_struts(w,l,h,strut_points,struts_jun_bi_pad)

%% identify all free skeletal brnahces

if size(strut_points,1)==1
    
    [a,b,c] = ind2sub([w,l,h],strut_points);
    a=a+1;b=b+1;c=c+1;
    window=struts_jun_bi_pad(a-1:a+1,b-1:b+1,c-1:c+1);
    if sum(window(:))==1 
        flag=1;
    else if sum(window(:))==2
            flag=2;
        else
            flag=0;
        end
    end

else
    [a1,b1,c1] = ind2sub([w,l,h],strut_points(1));
    a1=a1+1;b1=b1+1;c1=c1+1;
    [a2,b2,c2] = ind2sub([w,l,h],strut_points(end));
    a2=a2+1;b2=b2+1;c2=c2+1;
    window1=struts_jun_bi_pad(a1-1:a1+1,b1-1:b1+1,c1-1:c1+1);
    window2=struts_jun_bi_pad(a2-1:a2+1,b2-1:b2+1,c2-1:c2+1);
    if sum(window1(:))==1 && sum(window2(:))==0
        flag=1;
    else if sum(window1(:))==0 && sum(window2(:))==1
            flag=1;
    else if sum(window1(:))==1 && sum(window2(:))==1
            flag=2;
        else
            flag=0;
        end
        end
    end
    
    
end




end