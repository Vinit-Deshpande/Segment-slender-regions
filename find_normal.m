function normal_avg=find_normal(strut_points,x,y,z,j)

if length(strut_points)<4
    no=length(strut_points);
    else if length(strut_points)<10
         no=round(length(strut_points)/2);
        else if length(strut_points)<20
              no=round(length(strut_points)/3);
            else
             no=round(length(strut_points)/4);
            end
        end
end

S = interpc([x,y,z],no,'spline');
S_new = interpc([S(:,1),S(:,2),S(:,3)],100,'spline');



    jn = dsearchn(S_new,[x(j) y(j) z(j)]);
    if jn ==1
        id=[0 1 2 3 4 5 6];
    else if jn ==2
            id=[-1 0 1 2 3 4 5];
        else if jn ==3
                id=[-2 -1 0 1 2 3 4];
            else if jn >3 && jn <size(S_new,1)-2
                    id=[-3 -2 -1 0 1 2 3];
                else if jn ==size(S_new,1)-2
                        id=[-4 -3 -2 -1 0 1 2];
                    else if jn ==size(S_new,1)-1
                            id=[-5 -4 -3 -2 -1 0 1];
                        else if jn ==size(S_new,1)
                             id=[-6 -5 -4 -3 -2 -1 0];
                            end
                         end
                    end
                end
            end
        end
    end
    
         
          normal1=[S_new(jn+id(2),1)-S_new(jn+id(1),1),S_new(jn+id(2),2)-S_new(jn+id(1),2),S_new(jn+id(2),3)-S_new(jn+id(1),3)];
          normal1=normal1/(norm(normal1));
          normal2=[S_new(jn+id(3),1)-S_new(jn+id(2),1),S_new(jn+id(3),2)-S_new(jn+id(2),2),S_new(jn+id(3),3)-S_new(jn+id(2),3)];
          normal2=normal2/(norm(normal2));
          normal3=[S_new(jn+id(4),1)-S_new(jn+id(3),1),S_new(jn+id(4),2)-S_new(jn+id(3),2),S_new(jn+id(4),3)-S_new(jn+id(3),3)];
          normal3=normal3/(norm(normal3));
          normal4=[S_new(jn+id(5),1)-S_new(jn+id(4),1),S_new(jn+id(5),2)-S_new(jn+id(4),2),S_new(jn+id(5),3)-S_new(jn+id(4),3)];
          normal4=normal4/(norm(normal4));
          normal5=[S_new(jn+id(6),1)-S_new(jn+id(5),1),S_new(jn+id(6),2)-S_new(jn+id(5),2),S_new(jn+id(6),3)-S_new(jn+id(5),3)];
          normal5=normal5/(norm(normal5));
          normal6=[S_new(jn+id(7),1)-S_new(jn+id(6),1),S_new(jn+id(7),2)-S_new(jn+id(6),2),S_new(jn+id(7),3)-S_new(jn+id(6),3)];
          normal6=normal6/(norm(normal6));  
          
          normal=[normal1;normal2;normal3;normal4;normal5;normal6];
          normal_avg=mean(normal);
end


    
    
    
       
        