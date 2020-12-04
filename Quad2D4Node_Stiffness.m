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