function [ minH, maxH, minAng, maxAng, numVert ] = meshInfo( nodes, poly )
%meshInfo - gets information about the quality of the mesh
%
%inputs: nodes - node list
%poly - list of triangles in mesh
%
%Outputs: minH - minimum side length
%maxH - Max side length
%minAng - minimum angle in radians
%maxAng - maximum angle in radians
%numVert - number of vertices

numVert = size(nodes,1);

maxH = 0;
minH = Inf;
minAng = Inf;
maxAng = 0;

for i=1:size(poly,1)%for every polygon
   a = cell(3,1);
   for j=1:3%for every edge
      a{j} = nodes{poly(i,j),1} - nodes{poly(i,mod(j+1,3)+1),1};
      
      %update minH and maxH
      if(norm(a{j}) > maxH)
         maxH = norm(a{j});
      elseif(norm(a{j}) < minH)
         minH = norm(a{j});
      end
   end
   
   for j=1:3%for every angle
      %calculate angle
      temp = acosd(abs(sum(a{j}.*a{mod(j+1,3)+1}))/(norm(a{j})*norm(a{mod(j+1,3)+1})));
      
      %update minAng and maxAng
       %update minH and maxH
      if(temp > maxAng)
         maxAng = temp;
      elseif(temp < minAng)
         minAng = temp;
      end
   end
end

end

