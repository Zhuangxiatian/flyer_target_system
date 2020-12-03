clc
clear all

Ep='Youngs modulus is:';
E=input(Ep);
NUp='Poissons ratio is:';
NU=input(NUp);
hp='Thickness is:';
h=input(hp);
lep='Element size is:';
le=input(lep);
Nelxp='The number of element in x-direction is:';
Nelx=input(Nelxp);
Nelyp='The number of element in y-direction is:';
Nely=input(Nelyp);
rho2p='Density of flyer is:';
rho2=input(rho2p);
rho1p='Density of target is:';
rho1=input(rho1p);
Bop='The number of flyer element in x-direction:';
Bo=input(Bop);
dtp='The time step is:';
dt=input(dtp);
Tp='Final time is:';
T=input(Tp);
nfp='The number of Figure is:';
nf=input(nfp);
vinxp='The initial velocity in x-direction of flyer is:';
vinx=input(vinxp);
vinyp='The initial velocity in y-direction of flyer is:';
viny=input(vinyp);
dispxp='The initial displacement in x-direction of flyer is:';
dispx=input(dispxp);
dispyp='The initial displacement in y-direction of flyer is:';
dispy=input(dispyp);


%%
% set coordinate
coordinate=zeros((Nelx+1)*(Nely+1),2);
for y=1:(Nely+1)
    for x=1:(Nelx+1)
coordinate((y-1)*(Nelx+1)+x,1)=(x-1)*le;
coordinate((y-1)*(Nelx+1)+x,2)=(y-1)*le;
    end
end
%%
% set stiffness/mass matrix and damping ratio
KK=zeros(2*(Nelx+1)*(Nely+1),2*(Nelx+1)*(Nely+1));
MM=zeros(2*(Nelx+1)*(Nely+1),2*(Nelx+1)*(Nely+1));
[k1,k2,k3,massele1,massele2]=simplifiedprocess(E,NU,h,le,rho1,rho2);
for y=1:Nely
    for x=1:Nelx
        m=(Nelx+1)*(y-1)+x;
        p=(Nelx+1)*(y-1)+x+1;
        j=(Nelx+1)*y+x;
        i=(Nelx+1)*y+x+1;
        node((y-1)*Nelx+x,1)=m;
        node((y-1)*Nelx+x,2)=j;
        node((y-1)*Nelx+x,3)=i;
        node((y-1)*Nelx+x,4)=p;
        if x>Bo
            k=k2;
            massele=massele2;
        else
            k=k1;
            massele=massele1;
        end
        %z=Quad2D4Node_Assembly(KK,k,i,j,m,p);
        DOF(1)=2*i-1;
        DOF(2)=2*i;
        DOF(3)=2*j-1;
        DOF(4)=2*j;
        DOF(5)=2*m-1;
        DOF(6)=2*m;
        DOF(7)=2*p-1;
        DOF(8)=2*p;
for n1=1:8
    for n2=1:8
        KK(DOF(n1),DOF(n2))=KK(DOF(n1),DOF(n2))+k(n1,n2);
        MM(DOF(n1),DOF(n2))=MM(DOF(n1),DOF(n2))+massele(n1,n2);
    end
end    
    end
end

KK1=KK;

for y=1:Nely
    for x=1:Nelx
        m=(Nelx+1)*(y-1)+x;
        p=(Nelx+1)*(y-1)+x+1;
        j=(Nelx+1)*y+x;
        i=(Nelx+1)*y+x+1;
        if x>(Bo+1)
            k=k2;
            massele=massele2;
        else if x==(Bo+1)
            k=k3;
        else
            k=k1;
            massele=massele1;
            end
        end
        %z=Quad2D4Node_Assembly(KK,k,i,j,m,p);
        DOF(1)=2*i-1;
        DOF(2)=2*i;
        DOF(3)=2*j-1;
        DOF(4)=2*j;
        DOF(5)=2*m-1;
        DOF(6)=2*m;
        DOF(7)=2*p-1;
        DOF(8)=2*p;
for n1=1:8
    for n2=1:8
        KK(DOF(n1),DOF(n2))=KK(DOF(n1),DOF(n2))+k(n1,n2);
        MM(DOF(n1),DOF(n2))=MM(DOF(n1),DOF(n2))+massele(n1,n2);
    end
end
    end
end
KK2=KK;
C=zeros(2*(Nelx+1)*(Nely+1),2*(Nelx+1)*(Nely+1));
%%
% set initial condition of displacement, velocity and force
nt=T/dt;
D0=zeros(2*(Nelx+1)*(Nely+1),1);
V0=zeros(2*(Nelx+1)*(Nely+1),1);
F=zeros(2*(Nelx+1)*(Nely+1),nt+1);
for i=1:(Nely+1)
    V0(((i-1)*2*(Nelx+1)+2*(Bo+1)+1):2:((i-1)*2*(Nelx+1)+2*Nelx+1),1)=vinx;
    V0(((i-1)*2*(Nelx+1)+2*(Bo+1)+2):2:((i-1)*2*(Nelx+1)+2*(Nelx+1)),1)=viny;
end

for i=1:(Nely+1)
    D0(((i-1)*2*(Nelx+1)+2*(Bo+1)+1):2:((i-1)*2*(Nelx+1)+2*Nelx+1),1)=dispx;
    D0(((i-1)*2*(Nelx+1)+2*(Bo+1)+2):2:((i-1)*2*(Nelx+1)+2*(Nelx+1)),1)=dispy;
end
%%
% set boundary condition
bcdof(2:2:2*(Nelx+1))=1;
bcdof(2*Nely*(Nelx+1)+2:2:2*(Nely+1)*(Nelx+1))=1;
bcdof(2*Nely*(Nelx+1)+1:2*Nely*(Nelx+1)+2*(Bo+1))=1;
%%

kk=KK1;
cc=C;
mm=MM;
fd=F;
q0=D0;
dq0=V0;

[sdof,n2]=size(kk); 
 
dsp=zeros(sdof,nt);                                         % displacement matrix 
vel=zeros(sdof,nt);                                             % velocity matrix 
acc=zeros(sdof,nt);                                          % acceleration matrix 
 
dsp(:,1)=q0;                                                % initial displacement 
vel(:,1)=dq0;                                                   % initial velocity 
%-------------------------------------------------------------------------- 
%  (2) central difference scheme for time integration 
%-------------------------------------------------------------------------- 
acc(:,1)=(mm)\(fd(:,1)-kk*dsp(:,1)-cc*vel(:,1)); 
                                           % compute the initial acceleration (t=0) 
dsp0=dsp(:,1)-vel(:,1)*dt+0.5*acc(:,1)*dt^2; 
                                    % compute the fictitious displacement at time -dt 
ekk=mm/dt^2+cc/(2*dt); 
                                           % compute the effective stiffness matrix 
%------------------------------------------ 
%  (2.1) first step of the central difference scheme 
%------------------------------------------ 
efd=fd(:,1)-(kk-2*mm/dt^2)*dsp(:,1)-(mm/dt^2-cc/(2*dt))*dsp0; 
                                            %  compute the effective force vector 
dsp(:,1+1)=inv(ekk)*efd; 
                                                         % find the dsp at 1+dt 
for i=1:sdof                  % assign zero to acc, vel, dsp of the dofs associated with bc 
  if bcdof(i)==1 
    dsp(i,1)=0; 
    dsp(i,2)=0;
    vel(i,1)=0;
    acc(i,1)=0; 
  end 
end 
%------------------------------------------------- 
%  (2.2) subsequent steps of the central difference scheme 
%------------------------------------------------- 

count=zeros(nt);
for it=2:nt                                   % loop for each time step after first step 
  efd=fd(:,it)-(kk-2*mm/dt^2)*dsp(:,it)-(mm/dt^2-cc/(2*dt))*dsp(:,it-1); 
                                            %  compute the effective force vector 
  
      if vel(2*Bo+1)<vel(2*Bo+3)
          kk=KK2;
          count(it)=1;
      end
  
          
  dsp(:,it+1)=inv(ekk)*efd;                                   % find the dsp at t+dt 
  acc(:,it)=(dsp(:,it+1)-2*dsp(:,it)+dsp(:,it-1))/dt^2;                  % find the acc at t 
  vel(:,it)=(dsp(:,it+1)-dsp(:,it-1))/(2*dt);                           % find the vel at t 
  
  for i=1:sdof                % assign zero to acc, vel, dsp of the dofs associated with bc 
    if bcdof(i)==1 
      dsp(i,it)=0; 
      dsp(i,it+1)=0; 
      vel(i,it)=0; 
      acc(i,it)=0; 
    end 
  end 
 
end 
 
  acc(:,it+1)=acc(:,it); vel(:,it+1)=vel(:,it); 
 
if cc(1,1)==0 
  disp('The transient response results of undamping system') 
else 
  disp('The transient response results of damping system') 
end

vel=-transpose(vel);
dsp=transpose(dsp);
acc=transpose(acc);

gapf=round(nt/nf);

for i=1:gapf:nt
    for j=1:(Nelx+1)*(Nely+1)
    coordinatef(j,1)=coordinate(j,1)+dsp(i,(j-1)*2+1);
    coordinatef(j,2)=coordinate(j,2)+dsp(i,j*2);
    end
    velxf=vel(i,1:2:2*(Nelx+1)*(Nely+1)-1);
    PlotFieldonMesh(coordinatef,node,velxf);
end



function [acc,vel,dsp]=TransResp1(kk,cc,mm,fd,bcdof,nt,dt,q0,dq0,Nelx,Nely) 
%-------------------------------------------------------------------------- 
%  Purpose: 
%     The function subroutine TransResp1.m calculates transient response of 
%     a structural system using central difference scheme. 
%  Synopsis: 
%     [acc,vel,dsp]=TransResp1(kk,cc,mm,fd,bcdof,nt,dt,q0,dq0) 
%  Variable Description: 
%     Input parameters 
%       kk, cc, mm - stiffness, damping and mass matrices 
%       fd - Input or forcing influence matrix 
%       bcdof - Boundary condition dofs vector 
%       nt - Number of time steps 
%       dt - Time step size 
%       q0, dq0 - Initial condition vectors 
%     Output parameters 
%       acc - Acceleration response 
%       vel - Velocity response 
%       dsp - Displacement response 
%-------------------------------------------------------------------------- 
%  (1) initial condition 
%-------------------------------------------------------------------------- 
[sdof,n2]=size(kk); 
 
dsp=zeros(sdof,nt);                                         % displacement matrix 
vel=zeros(sdof,nt);                                             % velocity matrix 
acc=zeros(sdof,nt);                                          % acceleration matrix 
 
dsp(:,1)=q0;                                                % initial displacement 
vel(:,1)=dq0;                                                   % initial velocity 
%-------------------------------------------------------------------------- 
%  (2) central difference scheme for time integration 
%-------------------------------------------------------------------------- 
acc(:,1)=(mm)\(fd(:,1)-kk*dsp(:,1)-cc*vel(:,1)); 
                                           % compute the initial acceleration (t=0) 
dsp0=dsp(:,1)-vel(:,1)*dt+0.5*acc(:,1)*dt^2; 
                                    % compute the fictitious displacement at time -dt 
ekk=mm/dt^2+cc/(2*dt); 
                                           % compute the effective stiffness matrix 
%------------------------------------------ 
%  (2.1) first step of the central difference scheme 
%------------------------------------------ 
efd=fd(:,1)-(kk-2*mm/dt^2)*dsp(:,1)-(mm/dt^2-cc/(2*dt))*dsp0; 
                                            %  compute the effective force vector 
dsp(:,1+1)=inv(ekk)*efd; 
                                                         % find the dsp at 1+dt 
for i=1:sdof                  % assign zero to acc, vel, dsp of the dofs associated with bc 
  if bcdof(i)==1 
    dsp(i,1)=0; 
    dsp(i,2)=0;
    vel(i,1)=0;
    acc(i,1)=0; 
  end 
end 
%------------------------------------------------- 
%  (2.2) subsequent steps of the central difference scheme 
%------------------------------------------------- 
for it=2:nt                                   % loop for each time step after first step 
  efd=fd(:,it)-(kk-2*mm/dt^2)*dsp(:,it)-(mm/dt^2-cc/(2*dt))*dsp(:,it-1); 
                                            %  compute the effective force vector 
  %for xx=1:Nely+1
      %if vel(2*(xx-1)*(Nelx+1)+35)<vel(2*(xx-1)*(Nelx+1)+37)
          %efd(2*(xx-1)*(Nelx+1)+35)=0;
          %efd(2*(xx-1)*(Nelx+1)+37)=0;
      %end
  %end
          
  dsp(:,it+1)=inv(ekk)*efd;                                   % find the dsp at t+dt 
  acc(:,it)=(dsp(:,it+1)-2*dsp(:,it)+dsp(:,it-1))/dt^2;                  % find the acc at t 
  vel(:,it)=(dsp(:,it+1)-dsp(:,it-1))/(2*dt);                           % find the vel at t 
  
  for i=1:sdof                % assign zero to acc, vel, dsp of the dofs associated with bc 
    if bcdof(i)==1 
      dsp(i,it)=0; 
      dsp(i,it+1)=0; 
      vel(i,it)=0; 
      acc(i,it)=0; 
    end 
  end 
 
end 
 
  acc(:,it+1)=acc(:,it); vel(:,it+1)=vel(:,it); 
 
if cc(1,1)==0 
  disp('The transient response results of undamping system') 
else 
  disp('The transient response results of damping system') 
end
end


function k=Quad2D4Node_Stiffness(E,NU,h,xi,yi,xj,yj,xm,ym,xp,yp,ID)
 syms s t;
 a=(yi*(s-1)+yj*(-1-s)+ym*(1+s)+yp*(1-s))/4;
 b=(yi*(t-1)+yj*(1-t)+ym*(1+t)+yp*(-1-t))/4;
 c=(xi*(t-1)+xj*(1-t)+xm*(1+t)+xp*(-1-t))/4;
 d=(xi*(s-1)+xj*(-1-s)+xm*(1+s)+xp*(1-s))/4;
 B1=[a*(t-1)/4-b*(s-1)/4 0;0 c*(s-1)/4-d*(t-1)/4;c*(s-1)/4-d*(t-1)/4 a*(t-1)/4-b*(s-1)/4];%%
 B2=[a*(-t+1)/4-b*(-s-1)/4 0;0 c*(-s-1)/4-d*(-t+1)/4;c*(-s-1)/4-d*(-t+1)/4 a*(-t+1)/4-b*(-s-1)/4];%%
 B3=[a*(t+1)/4-b*(s+1)/4 0;0 c*(s+1)/4-d*(t+1)/4;c*(s+1)/4-d*(t+1)/4 a*(t+1)/4-b*(s+1)/4];%%
 B4=[a*(-t-1)/4-b*(-s+1)/4 0;0 c*(-s+1)/4-d*(-t-1)/4;c*(-s+1)/4-d*(-t-1)/4 a*(-t-1)/4-b*(-s+1)/4];%%
 Bfirst=[B1 B2 B3 B4];
 Jfirst=[0 1-t t-s s-1;t-1 0 s+1 -s-t;s-t -s-1 0 t+1;1-s s+t -t-1 0];
 J=[xi xj xm xp]*Jfirst*[yi;yj;ym;yp]/8;
 B=Bfirst/J;
 if ID==2
    D=(E/(1+NU)/(1-NU))*[1 NU 0;NU 1 0;0 0 (1-NU)/2];
 elseif ID==1
    D=1000000000*[32.699 1.257 0;1.257 4.265 0;0 0 1.432];
 end
 BD=J*transpose(B)*D*B;
 r=int(int(BD,t,-1,1),s,-1,1);
 z=h*r;%% h is the thickness
 k=double(z);
end


function massele=Quad2D4Node_Mass(rho1,rho2,h,le,ID)
if ID==1
    rho=rho1;
else
    rho=rho2;
end
    
massele=0.25*rho*h*le^2*[1 0 0 0 0 0 0 0;
                    0 1 0 0 0 0 0 0;
                    0 0 1 0 0 0 0 0;
                    0 0 0 1 0 0 0 0;
                    0 0 0 0 1 0 0 0;
                    0 0 0 0 0 1 0 0;
                    0 0 0 0 0 0 1 0;
                    0 0 0 0 0 0 0 1];
end

function [k1,k2,k3,massele1,massele2]=simplifiedprocess(E,NU,h,le,rho1,rho2)
x=1;
y=1;
xi=(x)*le;
yi=(y)*le;
xp=(x)*le;
yp=(y-1)*le;
xj=(x-1)*le;
yj=y*le;
xm=(x-1)*le;
ym=(y-1)*le;
ID=1;
k1=Quad2D4Node_Stiffness(E,NU,h,xi,yi,xj,yj,xm,ym,xp,yp,ID);
massele1=Quad2D4Node_Mass(rho1,rho2,h,le,ID);
ID=2;
k2=Quad2D4Node_Stiffness(E,NU,h,xi,yi,xj,yj,xm,ym,xp,yp,ID);
massele2=Quad2D4Node_Mass(rho1,rho2,h,le,ID);
k3=zeros(8,8);
end

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
% Dimensions of the colorbar     
% cpos = get(cbar,'position'); 
% cpos(3) = cpos(3)/4 ;   % Reduce the width of colorbar by half
% cpos(2) = cpos(2)+.1 ;
% cpos(4) = cpos(4)-.2 ;
%set(cbar,'Position',cpos) ;
brighten(0.5); 
     
% Title of the colorbar
set(get(cbar,'title'),'string','VAL');
%locate = get(cbar,'title');
%tpos = get(locate,'position');
%tpos(3) = tpos(3)+5. ;
%set(locate,'pos',tpos);

% Setting the values on colorbar
%
% get the color limits
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
