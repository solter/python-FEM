function [ nodes, newPoly ] = splitMesh( nodes, polys )
%splitMesh Takes a mesh and splits each edge into two
%
%inputs: nodes - n by 2 cell array of nodes. 1st cell is a length 2 array
%holding x and y values, second cell is a boolean true iff node is boundary
%
%polys - n by 3 array of integers listing the node number of each corner of
%vertex
%
%outputs:nodes - same as above spec but for split mesh
%newPoly - same as above spec but for split mesh

newNodes = cell(size(nodes,1));
newPoly = zeros(4*size(polys,1),3);

%store which nodes go to which polygons
for i=1:size(polys,1)
    
    %construct triangles which have 1 node from old triangle
    for j=1:3
        %the old node in the new triangle
        newPoly(4*i-j,j) = polys(i,j);
        
        %where the new nodes will be in the new triangle
        for k = 1:2
            if(polys(i,j) < polys(i,1+mod(j+k-1,3)))%make sure that no duplicates are created
              newNodes{polys(i,j),polys(i,1+mod(j+k-1,3))} = [newNodes{polys(i,j),polys(i,1+mod(j+k-1,3))}; [4*i-j,1+mod(j+k-1,3)]];
            else
              newNodes{polys(i,1+mod(j+k-1,3)),polys(i,j)} = [newNodes{polys(i,1+mod(j+k-1,3)),polys(i,j)}; [4*i-j,1+mod(j+k-1,3)]];
            end
        end
    end
    
    %%%%%%%%%%%%%%The above loop is functionally equivalent to the
    %%%%%%%%%%%%%%following except it makes sure that the first index in
    %%%%%%%%%%%%%%the nodes is lower than the second    
%     %construct first new triangle:
%     %the old node in the new triangle
%     newPoly(4*i-1,1) = polys(i,1);
%     %where the new nodes will be in the new triangle
%     newNodes{polys(i,1),polys(i,2)} = [newNodes{polys(i,1),polys(i,2)}; [4*i-1,2]];
%     newNodes{polys(i,1),polys(i,3)} = [newNodes{polys(i,1),polys(i,3)}; [4*i-1,3]];
%     
%     %construct second new triangle:
%     %the old node in the new triangle
%     newPoly(4*i-2,2) = polys(i,2);
%     %where the new nodes will be in the new triangle
%     newNodes{polys(i,2),polys(i,3)} = [newNodes{polys(i,2),polys(i,3)}; [4*i-2,3]];
%     newNodes{polys(i,2),polys(i,1)} = [newNodes{polys(i,2),polys(i,1)}; [4*i-2,1]];
%     
%     %construct second new triangle:
%     %the old node in the new triangle
%     newPoly(4*i-3,3) = polys(i,3);
%     %where the new nodes will be in the new triangle
%     newNodes{polys(i,3),polys(i,1)} = [newNodes{polys(i,3),polys(i,1)}; [4*i-3,1]];
%     newNodes{polys(i,3),polys(i,2)} = [newNodes{polys(i,3),polys(i,2)}; [4*i-3,2]];

    %construct center triangle which had no old nodes in it
    for j=1:3
        if(polys(i,j) < polys(i,1+mod(j+1-1,3)))
            newNodes{polys(i,j),polys(i,1+mod(j+1-1,3))} = [newNodes{polys(i,j),polys(i,1+mod(j+1-1,3))}; [4*i,j]];
        else
            newNodes{polys(i,1+mod(j+1-1,3)),polys(i,j)} = [newNodes{polys(i,1+mod(j+1-1,3)),polys(i,j)}; [4*i,j]];
        end
    end
end

%move all the nodes in the newNodes cell array into nodes
origSize = size(nodes,1);
for i=1:origSize
   for j=i:origSize
       if(isempty(newNodes{i,j})) continue; end%skip if there are no nodes at this connection
       
       pos = (nodes{i,1}+nodes{j,1})/2;
       if(size(newNodes{i,j},1) == 3)%if it is a physical boundary
         if(nodes{i,2}(1) == -1 || nodes{j,2}(1)==-1)
            bound = [-1 0];%use same boundary condition as adjacent node
         else
            bound = [1 0];
         end
       elseif(nodes{i,2}(1) == 2 && nodes{j,2}(1) == 2 && 2==1)%if it is an internally constant potential
          if(nodes{i,2}(2) == nodes{i,2}(2))%if same voltage
             bound = nodes{i,2};%make appropriate boundary
          else%if different voltage
             av = (nodes{i,2}(2) + nodes{j,2}(2))/2;%find average voltage
             %assign average voltage to all three nodes
             nodes{i,2}(2) = av;
             nodes{j,2}(2) = av;
             bound = [2, av];
          end
       else%if it is free to vary
          bound = [0 0];
       end
       nodes = vertcat(nodes,{pos bound});%add the midpoint to the node list
       
       n = size(nodes,1);
       
       for k=1:size(newNodes{i,j},1)%for every polygon each node is part of
          newPoly(newNodes{i,j}(k,1),newNodes{i,j}(k,2)) = n;%upddate polygon lists
       end
   end
end

end

