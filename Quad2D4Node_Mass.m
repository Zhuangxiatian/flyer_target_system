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