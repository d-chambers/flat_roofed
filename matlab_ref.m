

for j=1:20
    
    %Inputs
    
    T=1; % Beam thickness (m)
    S=j; % Span (m)
    E=1.0e10; % Intact rock Young's Modulus (Pa)
    dens=3000.0; % Intact rock density (kg/m^3)
    jknx=5.0e9; % Joint stiffness (Pa/m)
    sjx=1; % Spacing of vertical joints
    
    % Stiffness Calculations
    
    Ex=1./((1./E)+(1./(jknx.*sjx)));
    Erm=Ex;
    
    % Other Calculations
    
    g=(9.81*dens);
    
    Nmin=1;
    Nmax=0;
    N=0.01:0.01:1;
    N=N';
    Fmp=1E+100;
    mm=0;
    Np=0;
    disp=0;
    
    for i=1:100
        Zo=T*(1-(2/3)*N(i));
        L=S+(8*Zo*Zo/(3*S));
        del_l=0;
        del_l_prev=100;
        k=0;
        flag=0;
        while ((abs(del_l-del_l_prev))>0.000001)
            zchk=((8*Zo*Zo/(3*S))-del_l);
            if zchk>=0
                k=k+1;
                Z=sqrt(3*S*zchk/8);
                Fm=g*S*S/(4*N(i)*Z);
                Fav=Fm*((2/3)+N(i))/3;
                del_l_prev=del_l;
                del_l=(Fav/Erm)*L;
                mm=mm+1;
                
                bac(mm,1)=N(i);
                bac(mm,2)=Fm;
                bac(mm,3)=Z;
                bac(mm,4)=Zo;
            else
                del_l=del_l_prev;
                flag=1;
            end
            
            
        end
        
        if k>1
            if flag==0
                if N(i)<Nmin
                    Nmin=N(i);
                end
                if N(i)>Nmax
                    Nmax=N(i);
                end
            end
        end
        
        if Fmp>Fm
            if k>1
                if flag==0
                    
                    Fmp=Fm;
                    Np=N(i);
                    Zp=Z;
                    Zop=Zo;
                    disp=Zop-Zp;
                    
                end
            end
        end
        
        
        
    end
    
    % Outputs
    
    str(j,1)=i/100; % Number of iterations
    str(j,2)=Np; % N
    str(j,3)=100*(1-(Nmax-Nmin)); % Buckling Limit
    str(j,4)=-disp; % Midspan displacement
    str(j,5)=Fmp; % Fm
    str(j,6)=j; % Span
    
end