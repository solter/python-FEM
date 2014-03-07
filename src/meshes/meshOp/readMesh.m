function [ nodes , polys] = readMesh( filename , varargin)
%readMesh Takes in the string filename, and will import both the
%<filename>.ele and <filename>.node. Please Note: the files may not contain
%any numbers within their comments
%
%inputs: filename - string representing the file name
%2nd input - if present, boundaries will be assigned this value which is
%code for a boundary condition type instead of 1.
%
%Codes: -1,0 - neuman; 1 - dirichlet
%
%outputs: nodes - the nodes
%polys - the polygon mesh

if(length(varargin)==1)
    boundtype = varargin{1};
else
    boundtype = true;
end

%import the node file
fileID = fopen(strcat(filename,'.node'),'r','native');
nodeRaw = fscanf(fileID,'%f',Inf);
fclose(fileID);
valToSkip = nodeRaw(3) + nodeRaw(4);
if(nodeRaw(4) == 1)%if file indicates boundary values
boundVal =  true;
else
boundVal = false;
end

%find the number of nodes and the vertex count
numVert = nodeRaw(1,1);
nodes = cell(numVert,2);
for i = 1:numVert
	nodes{i,1} = [-1,-1];
end

i = 5;
while i < (size(nodeRaw,1))

    %if the node count is greater than number of vertices, throw error
    if(nodeRaw(i) > numVert)
        str = sprintf('%s.node has a node number which is \ngreater than the total number of nodes.',filename);
        throw(MException('FileIO:ImproperFormat',str));
    else%otherwise, copy over to sanitized node structure
        nodes{nodeRaw(i),1}(1) = nodeRaw(i+1);
        nodes{nodeRaw(i),1}(2) = nodeRaw(i+2);
        if(boundVal)
            if(nodeRaw(i + 2 + nodeRaw(3) + 1) == 1)
               nodes{nodeRaw(i),2} = [boundtype 0];
            else
               nodes{nodeRaw(i),2} =[false 0];
            end
        end
    end  
    
    i = i + 2 + valToSkip + 1;
 end

%if not all the nodes were defined
%if(min(nodes(:,1)) < 0)
 %    str = sprintf('%s.node has fewer nodes than it indicates it should have.',filename);
%    throw(MException('FileIO:ImproperFormat',str));
%end

%import the element file
fileID = fopen(strcat(filename,'.ele'),'r','native');
polyRaw = fscanf(fileID,'%f',Inf);
fclose(fileID);
valToSkip = polyRaw(3);
vc = polyRaw(2);%number of vertices per polygon

%find the number of nodes and the vertex count
numPoly = polyRaw(1,1);
polys = ones(numPoly,vc).*-1;

i = 4;
while i < (size(polyRaw,1))

    %if the node count is greater than number of vertices, throw error
    if(polyRaw(i) > numPoly)
        str = sprintf('%s.ele has a polygon number which is \ngreater than the total number of polygons.',filename);
        throw(MException('FileIO:ImproperFormat',str));
    else%otherwise, copy over to sanitized node structure
       for j = 1:vc
          polys(polyRaw(i), j) = polyRaw(i+j);
       end
    end  
    
    i = i + vc + valToSkip + 1;
end

%if not all the polygons were defined
if(min(polys(:,1)) < 0)
     str = sprintf('%s.ele has fewer polygons than it indicates it should have.',filename);
     throw(MException('FileIO:ImproperFormat',str));
end

end

