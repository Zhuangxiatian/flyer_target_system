function PlotFieldonMesh(coordinate,nodes,component)

nel = length(nodes) ;                  % number of elements
nnel = size(nodes,2);                % number of nodes per element
% 
% Initialization of the required matrices
X = zeros(nnel,nel) ;
Y = zeros(nnel,nel) ;
profile = zeros(nnel,nel) ;
%
for iel=1:nel   
     for i=1:nnel
     nd(i)=nodes(iel,i);         % extract connected node for (iel)-th element
     X(i,iel)=coordinate(nd(i),1);    % extract x value of the node
     Y(i,iel)=coordinate(nd(i),2);    % extract y value of the node
     end   
     profile(:,iel) = component(nd') ;         % extract component value of the node 
end
    
% Plotting the FEM mesh and profile of the given component
     f3 = figure ;
     set(f3,'name','Postprocessing','numbertitle','off') ;
     plot(X,Y,'k')
     fill(X,Y,profile)
     axis off ;
     % Colorbar Setting
     SetColorbar
end
 
function SetColorbar
cbar = colorbar;
brighten(0.5); 
set(get(cbar,'title'),'string','VAL');
clim = caxis;
ylim(cbar,[clim(1) clim(2)]);
numpts = 24 ;    % Number of points to be displayed on colorbar
kssv = linspace(clim(1),clim(2),numpts);
set(cbar,'YtickMode','manual','YTick',kssv); % Set the tickmode to manual
for i = 1:numpts
    imep = num2str(kssv(i),'%+3.2E');
    vasu(i) = {imep} ;
end
set(cbar,'YTickLabel',vasu(1:numpts),'fontsize',9);
end