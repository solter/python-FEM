function [] = writeMesh( nodes, poly, filename )
%writeMesh saves the mesh defined by the inputs into filename
%
%inputs: nodes - the vertices of the mesh
%poly - the polygons of the mesh
%filename - a string which will be saved into mesh format as filename.node
%and filename.ele

%save the nodes
bound = 0;
if(~isempty(nodes{1,2}))
	bound = 1;
end

fileID = fopen(strcat(filename,'.node'),'w');
fprintf(fileID,'%d\t2\t0\t%d\n',size(nodes,1),bound);

if(bound == 0)
for i = 1:size(nodes,1)
   fprintf(fileID,'%d\t%f\t%f\n',i,nodes{i,1}(1),nodes{i,1}(2));
end
else
for i = 1:size(nodes,1)
b = 0;
if(nodes{i,2}(1))
b = 1;
end
   fprintf(fileID,'%d\t%f\t%f\t%d\n',i,nodes{i,1}(1),nodes{i,1}(2),b);
end
	
end

%close the .node file
fclose(fileID);

%save data for polygons
fileID = fopen(strcat(filename,'.ele'),'w');
fprintf(fileID,'%d\t%d\t0',size(poly,1),size(poly,2));

for i = 1:size(poly,1)
   fprintf(fileID,'\n%d',i);
   for j = 1:size(poly,2)
       fprintf(fileID,'\t%d',poly(i,j));
   end
end

%construct and save the .ele file
fclose(fileID);

end

