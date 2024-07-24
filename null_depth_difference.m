clc;
clear;
close all;

%global cartesian co-ordinates
x_h=[1 0 0];
y_h=[0 1 0];
z_h=[0 0 1];



D=18;
D0=17;
f_D=0.4;
f=f_D*D;
v=1.5*10^9;

%Feed co-ordinate system
xfeed_h=x_h;
yfeed_h=-y_h;
zfeed_h=-z_h;

phi_p=0; %plane angle
theta_0=2*atan(1/(4*f_D));

%differenial surface (gridding dimensions)
drho=0.0207;
dphi=0.02;
drho1=0.1;
dphi1=0.1;

c=3e8;
eta=376.73;
i_count=1;
h_count=1;
j_count=1;
k_count=1;
q=1.14;
E_t=0;
theta_can=1.75; %nulling direction in deg
lambda=c/v;
beta=2*pi/lambda;
r_h0=[sind(theta_can)*cosd(phi_p) sind(theta_can)*sind(phi_p) cosd(theta_can)];

%Calculating the electric field to be nulled at 1.75 deg
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

        E_s = (J_s*exp(1i*beta*dot(r_h0,rfeed_v))*ds);
        E_t=E_t+E_s;
      
    end
    
end

%Calculating the phase value of cn to nullify the electric field obtained

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
    
        E_s1 = (J_s1.*exp(1i*beta*dot(r_h0,rfeed_v)).*ds);
        E_s2 = (J_s2.*exp(1i*beta*dot(r_h0,rfeed_v)).*ds);
        
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

%%Obtain the co-pol component for each frequency at 1.75 deg when nulling is implemented
for v=1400000000:100000:1600000000
    beta=2*pi*v/c;

    theta_1=1;
    r_h=[sind(theta_1)*cosd(phi_p) sind(theta_1)*sind(phi_p) cosd(theta_1)];
    r_mag = norm(r_h);
    E_t=[0 0 0];
    E_t0=[0 0 0];
    ecop_h=((-((1-cosd(theta_1))*sind(phi_p)*cosd(phi_p))).*x_h)+((1-sind(phi_p)*sind(phi_p)*(1-cosd(theta_1))).*y_h)-((sind(theta_1)*sind(phi_p)).*z_h);
    ecross_h=((1-cosd(phi_p)*cosd(phi_p)*(1-cosd(theta_1))).*x_h)-(((1-cosd(theta_1))*cosd(phi_p)*sind(phi_p)).*y_h)-(sind(theta_1)*cosd(phi_p).*z_h);
    
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
    
            E_s = (J_s*exp(1i*beta*dot(r_h0,rfeed_v))*ds);
            E_t=E_t+E_s;
          
        end
        
    end

    i_count=1;
    for rho=D0/2+drho:drho1:D/2+drho

        for phi=dphi:dphi1:2*pi+dphi

            phi_f=-phi;
            thetaf = -2*atan(rho/(2*f));
            rf=f*(sec(thetaf/2))^2;
            zf=-rf*cos(thetaf); %projection of the element on the +z axis
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
    E_t2=0;
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
    
            E_s = (J_s*exp(1i*beta*dot(r_h0,rfeed_v))*ds);
            E_t2=E_t2+E_s;
          
        end
    
    end

    E_dt=dot(ecop_h,E_t);
    E_dt2=dot(ecop_h,E_t2);
    E_xt=dot(ecross_h,E_t);
    E_xt2=dot(ecross_h,E_t2);
    copol_pattern0(h_count,1)=(4*norm(E_xt).^2)/0.0764; %Co-pol component in quiescent state
    copol_pattern2(h_count,1)=(4*norm(E_xt2).^2)/0.0764; %Co-pol component in ERS 1 bit phase

    h_count=h_count+1;
 
end

angle_value=1400000000:100000:1600000000;
copol=10*log10(copol_pattern2)-10*log10(copol_pattern0); %Difference in the depth of the null
plot(angle_value,copol,'b-')
grid on
xlabel('frequency of operation [Hz]')
ylabel('Diff in magnitude of co-pol pattern [dBi]')


