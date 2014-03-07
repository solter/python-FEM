function [  ] = plotObj( nodes, fig )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

figure(fig);
hold on;
for i=1:size(nodes,1)
    
    if(nodes{i,2}(1) == 2)
     plot(nodes{i,1}(1),nodes{i,1}(2),'k.');
    end
    
end

set(gca(),'dataaspectratio',[1,1,1])

end

