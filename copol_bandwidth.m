clc;
clear;
close all;

%Global co-ordinate system
x_h=[1 0 0];
y_h=[0 1 0];
z_h=[0 0 1];

D=18;
D0=17;
f_D=0.4;
f=f_D*D;
mu=4*pi*10^-7;

%Feed co-ordinate system
xfeed_h=x_h;
yfeed_h=-y_h;
zfeed_h=-z_h;

phi_p=0; %Plane angle
theta_0=2*atan(1/(4*f_D));
eta=376.73;
i_count=1;
h_count=1;
k_count=1;
m_count=1;
q=1.14;


for theta_can=1:0.01:2.8 %Range of the angle of observation
    theta_1=theta_can;
    r_h=[sind(theta_1)*cosd(phi_p) sind(theta_1)*sind(phi_p) cosd(theta_1)];
    ecop_h=((-((1-cosd(theta_1))*sind(phi_p)*cosd(phi_p))).*x_h)+((1-sind(phi_p)*sind(phi_p)*(1-cosd(theta_1))).*y_h)-((sind(theta_1)*sind(phi_p)).*z_h);
    v=1.5*10^9; %pre-defined center frequency
    c=3e8;  
    lambda=c/v;

    %differential surface area (gridding size)
    drho=0.1036*lambda;
    dphi=0.1*lambda;
    drho1=0.5*lambda;
    dphi1=0.5*lambda;

    beta=2*pi/lambda;
    k_count=1;
    j_count=1;
    i_count=1;
    E_t=0;

    %%Calculate the electric field which is to be nulled at the range of angle of observation
    for rho=drho:drho:D0/2+drho
           
        for phi=dphi:dphi:2*pi+dphi
            
            phi_f=-phi;
            thetaf = -2*atan(rho/(2*f));

            rf=f*(sec(thetaf/2))^2;
            zf=-rf*cos(thetaf); %projection of the element on the +z axis
            rfeed_v = [rho*cos(phi_f) rho*sin(phi_f) zf];
            rf=norm(rfeed_v);
          
            rf_h=rfeed_v/rf;
            E_i = (cross(cross(yfeed_h,rf_h),rf_h))*((exp(-1i*beta*rf))/rf);

            H_i = ((cross(rf_h,E_i)))*cos(thetaf)^q;
            
            ds = ((((4*f*f)+(rho*rho))^0.5)*rho*drho*dphi*cos(thetaf/2))/(2*f);             
            rho_h  = xfeed_h.*cos(phi_f)+yfeed_h.*sin(phi_f);
            n_h = ((-rho.*rho_h) +(2*f.*z_h))/(((4*f*f) + (rho*rho))^0.5);
            
            J_s = cross(2*n_h,H_i);

            E_s = (J_s*exp(1i*beta*dot(r_h,rfeed_v))*ds);
            E_t=E_t+E_s;
         
        end
        
    end
    
    E_dt=dot(ecop_h,E_t);
    E_t=E_dt;
    i_count=1;

    %%Obtain the value of cn which can reduce the electric field (Serial Search Algorithm)
    for rho=D0/2+drho:drho1:D/2+drho
        for phi=dphi:dphi1:2*pi+dphi
    
            phi_f=-phi;
            thetaf = -2*atan(rho/(2*f));
            rf=f*(sec(thetaf/2))^2;
            zf=-rf*cos(thetaf);
            rfeed_v = [rho*cos(phi_f) rho*sin(phi_f) zf];
            rf=norm(rfeed_v);
    
            rf_h=rfeed_v/rf;
            E_i = (cross(cross(yfeed_h,rf_h),rf_h))*((exp(-1i*beta*rf))/rf);
    
            H_i = ((cross(rf_h,E_i)))*cos(thetaf)^q;
    
            ds = ((((4*f*f)+(rho*rho))^0.5)*rho*drho1*dphi1)/(2*f);           
    
            rho_h  = xfeed_h.*cos(phi_f)+yfeed_h.*sin(phi_f);
            n_h = ((-rho.*rho_h) +(2*f.*z_h))/(((4*f*f) + (rho*rho))^0.5);
    
            cn1=1;
            cn2=-1;
    
            J_s = cross(2*n_h,H_i);
    
            J_s1=cn1.*J_s;
            J_s2=cn2.*J_s;
    
            E_s1 = (J_s1.*exp(1i*beta*dot(r_h,rfeed_v)).*ds);
            E_s2 = (J_s2.*exp(1i*beta*dot(r_h,rfeed_v)).*ds);
    
            E_t1=E_t+E_s1;
            E_t2=E_t+E_s2;
    
            if abs(E_t1(1,2))<abs(E_t2(1,2))
                cn=cn1;
                E_t=E_t1;
            else
                cn=cn2;
                E_t=E_t2;
            end
            cn_matrix(i_count,1)=cn;
            i_count=i_count+1;
        end
        
    end
    h_count=1;

    %%Calculate the co-pol component for a range of frequency to determine the bandwidth of the null
    for v=1400000000:1000000:1600000000
        lambda=c/v;
        beta=2*pi/lambda; 
        E_t=[0 0 0];
        E_t0=[0 0 0];
        
        for rho=drho:drho:D/2+drho
    
            for phi=dphi:dphi:2*pi+dphi
    
                phi_f=-phi;
                thetaf = -2*atan(rho/(2*f));
    
                rf=f*(sec(thetaf/2))^2;
                zf=-rf*cos(thetaf);
                rfeed_v = [rho*cos(phi_f) rho*sin(phi_f) zf];
                rf=norm(rfeed_v);
 
                rf_h=rfeed_v/rf;
                E_i = (cross(cross(yfeed_h,rf_h),rf_h))*((exp(-1i*beta*rf))/rf);
    
                H_i = ((cross(rf_h,E_i)))*cos(thetaf)^q;
    
                ds = ((((4*f*f)+(rho*rho))^0.5)*rho*drho*dphi*cos(thetaf/2))/(2*f);           
                rho_h  = xfeed_h.*cos(phi_f)+yfeed_h.*sin(phi_f);
                n_h = ((-rho.*rho_h) +(2*f.*z_h))/(((4*f*f) + (rho*rho))^0.5);
    
                J_s = cross(2*n_h,H_i);
    
                E_s = (J_s*exp(1i*beta*dot(r_h,rfeed_v))*ds);
                E_t0=E_t0+E_s;
    
            end
    
        end

        E_dt0=dot(ecop_h,E_t0);
        copol_pattern0(k_count,1)=(4*norm(E_dt0).^2)/0.0764;
        k_count=k_count+1;
    
        for rho=drho:drho:D0/2+drho
            
            for phi=dphi:dphi:2*pi+dphi
                
                phi_f=-phi;
                thetaf = -2*atan(rho/(2*f));
    
                rf=f*(sec(thetaf/2))^2;
                zf=-rf*cos(thetaf);
                rfeed_v = [rho*cos(phi_f) rho*sin(phi_f) zf];
                rf=norm(rfeed_v);
    
                
                rf_h=rfeed_v/rf;
                E_i = (cross(cross(yfeed_h,rf_h),rf_h))*((exp(-1i*beta*rf))/rf);
    
                H_i = ((cross(rf_h,E_i)))*cos(thetaf)^q;
                
                ds = ((((4*f*f)+(rho*rho))^0.5)*rho*drho*dphi*cos(thetaf/2))/(2*f);           
                rho_h  = xfeed_h.*cos(phi_f)+yfeed_h.*sin(phi_f);
                n_h = ((-rho.*rho_h) +(2*f.*z_h))/(((4*f*f) + (rho*rho))^0.5);
    
                J_s = cross(2*n_h,H_i);
    
                E_s = (J_s*exp(1i*beta*dot(r_h,rfeed_v))*ds);
                E_t=E_t+E_s;
              
            end
            
        end
    
    
        i_count=1;
        
        for rho=D0/2+drho:drho1:D/2+drho
    
            for phi=dphi:dphi1:2*pi+dphi
    
                phi_f=-phi;
                thetaf = -2*atan(rho/(2*f));
                rf=f*(sec(thetaf/2))^2;
                zf=-rf*cos(thetaf);
                rfeed_v = [rho*cos(phi_f) rho*sin(phi_f) zf];
                rf=norm(rfeed_v);
                
                rf_h=rfeed_v/rf;
                E_i = (cross(cross(yfeed_h,rf_h),rf_h))*((exp(-1i*beta*rf))/rf);
    
                H_i = ((cross(rf_h,E_i)))*cos(thetaf)^q;
                
                ds = ((((4*f*f)+(rho*rho))^0.5)*rho*drho1*dphi1)/(2*f);           
          
                rho_h  = xfeed_h.*cos(phi_f)+yfeed_h.*sin(phi_f);
                n_h = ((-rho.*rho_h) +(2*f.*z_h))/(((4*f*f) + (rho*rho))^0.5);
    
                J_s = cross(2*n_h,H_i);
    
                J_s1=cn_matrix(i_count,1).*J_s;
    
                E_s1 = (J_s1.*exp(1i*beta*dot(r_h,rfeed_v)).*ds);
    
                E_t=E_t+E_s1;
                i_count=i_count+1;
            end
        end
    
        E_dt=dot(ecop_h,E_t);
        copol_pattern(h_count,1)=(4*norm(E_dt).^2)/0.0764;
        h_count=h_count+1;

    end

    copol_diff=10*log10(copol_pattern0)-10*log10(copol_pattern); %Difference in the magnitude of co-pol in the angle of observation to 
    g=1;
    freq=linspace(1400000000,16000000000,201);
    for j_count=1:h_count-1
        if copol_diff(j_count,1)>30 %Choosing the values exceeding 30 dBi
            copol(g,1)=copol_diff(j_count,1);
            firstvalue(g,1)=freq(1,j_count);
            g=g+1;
        end
    end

    if g==1
        bw=0;
    else
        bw=firstvalue(g-1,1)-firstvalue(1,1); %Subtracting the first and last value of corresponding co-pol magnitude exceeding 30 dBi
        bandwidth(m_count)=bw;
        m_count=m_count+1;
    end
end

angle_value=linspace(1,2.8,146);
plot(angle_value,bandwidth)
grid on 
xlabel('direction of null [degrees]')
ylabel('bandwidth [Hz]')