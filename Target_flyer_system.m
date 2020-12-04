clc
clear all

Ep='Youngs modulus of flyer is:';
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
