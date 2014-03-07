function [] = plotMesh( nodes, poly, varargin )
%plotMesh Plots a mesh of the given nodes with the given polygons.
%
%inputs: nodes - list of nodes
%poly - list of polygons
%3rd argument - name of the mesh
%4th argument - if present, this should be a figure handle. This will be
%the figure on which these are plotted
%5th argument - if present, should be a boolean. if true, will plot full
%picture, if false will not plot boundary points or lines. If the 4th
%argument is 0 then it will generate a new plot

plt = true;
name = '';

if(length(varargin) == 1)
    name = varargin{1};
    figure();
elseif(length(varargin) == 2)
    figure(varargin{2})
elseif(length(varargin) == 3)
    if(varargin{2} == 0)
        figure()
    else
        figure(varargin{2})
    end
    
    plt = varargin{3};
else
    figure()
end

hold on;

%name plot and put points on graph
title(name);
for i=1:size(nodes,1)
    
    str = '';
    if(nodes{i,2}(1) == 0)
       continue;
       %str = '*';
    else
        if(nodes{i,2}(1) == -1)
            str = 'o';
        elseif(nodes{i,2}(1) == 1)
           str = 's';
        else
            continue;
        end
    end
    plot(nodes{i,1}(1),nodes{i,1}(2),str);
end

if(plt)
for i = 1:size(poly,1)%for each polygon
    for j = 1:size(poly,2)%for each side
        temp = zeros(2,2);
        if(j == 1)%if first node of polygon, draw first point connected to last point
            temp = [nodes{poly(i,size(poly,2)),1};nodes{poly(i,j),1}];
        else%otherwise, connnect it to the last point
            temp = [nodes{poly(i,j-1),1};nodes{poly(i,j),1}];
        end
        line(temp(:,1),temp(:,2));
    end
end
end

set(gca(),'dataaspectratio',[1,1,1])
title('circ - neumann, square - dirichlet');

end

