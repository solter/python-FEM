function [] = plotColors(triags,vert,u,varargin)
%@author Peter Solfest
%plotColors plots the solution on the domain as a 2-D grid with a color spectrum
%ranging from red to blue from max to least
%
%inputs: triags - list of triangles
%vert - list of vertices
%u - solution at each vertex
%Note: these inputs must come from the readMesh_oct and fem2Ddiri functions

   if(length(varargin) ==1)
      figure(varargin{1});
   else
      figure();
   end
	maxu = max(u);
	minu = min(u);
	%for each triangle
	for i = 1:length(triags)
		%save vertices to a usable form
		verts(1,:) = vert{triags(i,1),1}(:);
		verts(2,:) = vert{triags(i,2),1}(:);
		verts(3,:) = vert{triags(i,3),1}(:);
		fac = [1:3];
		
		%find the average value of u across the triangle
		meanu = (u(triags(i,1))+ u(triags(i,2))+ u(triags(i,3)))/3;
		%find how large the mean u is relative to the range of u
		if(maxu-minu > 0)
         uporp = (meanu - minu)/(maxu-minu);
		else
			uporp = 0;
         meanu = 0;
		end
		%draw the triangle in the correct color
      green = 0;%1/(1+abs(meanu));
		patch('Faces',fac,'Vertices',verts,'FaceColor',[.8*uporp+.2,green,1-uporp], 'EdgeColor', 'none');
	end
	
	%label the graph
	title(sprintf('Red = %f       -->      Blue = %f;',maxu, minu),'FontSize',20);
	%xlabel('x axis','FontSize',20);
	%ylabel('y axis','FontSize',20);
end