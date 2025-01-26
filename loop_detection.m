function [end_point_dist,arc_length]=loop_detection(points)

arc_length=0;
for i=1:size(points,1)-1
    
    arc_length=arc_length+(sqrt(sum((points(i+1,:)-points(i,:)).^2)));    
       
end

end_point_dist=(sqrt(sum((points(end,:)-points(1,:)).^2)));


end