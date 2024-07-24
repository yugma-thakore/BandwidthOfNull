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
v=1.5*10^9;
c=3e8;

%Feed co-ordinate system
xfeed_h=x_h;
yfeed_h=-y_h;
zfeed_h=-z_h;
phi_p=0; %Plane angle
theta_0=2*atan(1/(4*f_D));
lambda=c/v;

%differential surface area
drho=0.1036*lambda;
dphi=0.1*lambda;

beta=2*pi/lambda;
eta=376.73;
i_count=1;
h_count=1;
q=1.14;

for theta_1=0:0.002:2 %Angle of observation in degrees
    r_h=[sind(theta_1)*cosd(phi_p) sind(theta_1)*sind(phi_p) cosd(theta_1)];
    r_mag = norm(r_h);
    E_t=[0 0 0];
    ecop_h=((-((1-cosd(theta_1))*sind(phi_p)*cosd(phi_p))).*x_h)+((1-sind(phi_p)*sind(phi_p)*(1-cosd(theta_1))).*y_h)-((sind(theta_1)*sind(phi_p)).*z_h);
    ecross_h=((1-cosd(phi_p)*cosd(phi_p)*(1-cosd(theta_1))).*x_h)-(((1-cosd(theta_1))*cosd(phi_p)*sind(phi_p)).*y_h)-(sind(theta_1)*cosd(phi_p).*z_h);

    for rho=drho:drho:D0/2+drho

        for phi=dphi:dphi:2*pi

            phi_f=-phi;
            thetaf = -2*atan(rho/(2*f));
            thetaf_h=xfeed_h*cos(thetaf)*cos(phi_f)+yfeed_h*cos(thetaf)*sin(phi_f)-zfeed_h*sin(thetaf);
            rf=f*(sec(thetaf/2))^2;
            zf=-rf*cos(thetaf);
            rfeed_v = [rho*cos(phi_f) rho*sin(phi_f) zf];
            rf=norm(rfeed_v);


            rf_h=rfeed_v/rf;
            E_i = (cross(cross(yfeed_h,rf_h),rf_h))*((exp(-1i*beta*rf))/rf);

            H_i = ((cross(rf_h,E_i)))*cos(thetaf)^q;

            ds = ((((4*f*f)+(rho*rho))^0.5)*rho*drho*cos(thetaf/2)*dphi)/(2*f);           

            rho_h  = xfeed_h.*cos(phi_f)+yfeed_h.*sin(phi_f);
            n_h = ((-rho.*rho_h) +(2*f.*z_h))/(((4*f*f) + (rho*rho))^0.5);

            J_s = cross(2*n_h,H_i);

            E_s = (J_s*exp(1i*beta*dot(r_h,rfeed_v))*ds);
            E_t=E_t+E_s;

        end

    end
    E_dt=dot(ecop_h,E_t);
    copol_pattern(i_count,1)=(norm(E_dt).^2)/0.0191;
    
    i_count=i_count+1;

end

angle_value=0:0.002:2;
plot(angle_value,10*log10(copol_pattern),'k-')
grid on

xlabel('Angle from reflector axis of rotation [deg]')
ylabel('Pattern [dBi]')
legend('D=18m')