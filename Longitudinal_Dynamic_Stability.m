%longitudinal Dynamic Stability  of Airplane E (case 1- power approach)

clear clc;
close all;
clear all;

%Data 

    %Flight Condition
alt = 0;%feet
a_dens = 0.002378; %slug/feet3
speed = 230;%fps
gravity_center = 0.29; % x_cg
attitude_ini = 11.7;% in rad
g_0 = 9.81;% m/s2

    %Geometry and Inertias
wing_area = 530; %feet2
wing_span = 38.7;%feet
wing_mean_chord = 16;% feet
weight = 33200;%lbs
i_xx = 23700;%slug feet2
i_yy = 117500;%slug feet2
i_zz = 133700;%slug feet2
i_xz = 1600; %slug feet2
    %Steady State Coefficients
c_L = 1;
c_D = 0.2;
c_T_X = 0.2;
c_m_ = 0;
c_m_T = 0;
    %longitudinal-Directional Derivatives
%Longitudinal Derivatives
Cmu = 0;
Cmalpha = -0.098;
Cmalphapoint = -0.95;
Cmq = -2;
CmTu = 0;
CmTalpha = 0;
CLu = 0;
CLalpha = 2.8;
Clalphapoint = 0;
CLq = 0;
CDalpha = 0.555;
CDu = 0;
CTXu = 0;
CLiH = -.024;
CDiH = -0.14;
CmiH = -0.322;

CLdeltaT = 0 ;
CDdeltaT = 0 ;
CMdeltaT = 0 ;

ft = 0.3048;
lb = 4.448; %kg*gravity
slug = 14.62;
ft2 = 0.0929;
slugft2 = 1.3558;
slugft3 = 515.383;

%Data conversion
    %Flight Condition
alt = alt*ft;%m
a_dens = a_dens*slugft3; %kg/m3
speed = speed*ft;%m/s
M = speed/340.3;

    %Geometry and Inertias
wing_area = wing_area*(ft2); %m2
wing_span = wing_span*ft; %m
wing_mean_chord = wing_mean_chord*ft; %m
weight = weight*lb/g_0; %kg
i_xx = i_xx*slugft2; %kg m2
i_yy = i_yy*slugft2; %kg m2
i_zz = i_zz*slugft2; %kg m2
i_xz = i_xz*slugft2; %kg m2
A=wing_span^2/wing_area;

%1.	Equation of longitudinal motion.
        %X=A*x+B*eta;
disp('1.	Equations of longitudinal motion.');
disp('X=A*x+B*eta');
% Calculation of A Matrix composant

 %flight condition
q=(1/2)*a_dens*speed^2;
coeff=(q*wing_area)/(weight*speed);
coeff2=(q*wing_area*wing_mean_chord)/(speed);

X_u=-coeff*(2*c_D+M*CDu);
X_w=coeff*(c_L-((2*c_L*CLalpha)/(pi*0.84*A)));
Z_u=-coeff*(2*c_L+(M^2*c_L)/(1-M^2));
Z_w=-coeff*(c_D+CLalpha);
Z_w_point = ((q*wing_area*wing_mean_chord)/(2*weight*(speed^2)))*(c_D+Clalphapoint);
Z_q=coeff2*(CLq)/(2*weight);
M_u=coeff2*M*Cmu/i_yy;
M_w=coeff2*Cmalpha/i_yy;
M_w_point=coeff2*(wing_mean_chord*Cmalphapoint)/(2*i_yy*speed);
M_q=coeff2*(wing_mean_chord*Cmq)/(2*i_yy);


%2.	The matrix A of aircraft

A=[X_u,X_w,0,-g_0*cos(attitude_ini);
   Z_u,Z_w,speed,-g_0*sin(attitude_ini);
   M_u+(M_w_point*Z_u),M_w+(M_w_point*Z_w),M_q+(speed*M_w_point),-M_w_point*g_0*sin(attitude_ini);
   0,0,1,0];
disp('A = ')
disp(A)

% Calculation of B Matrix composant

Xdeltae=coeff*CDiH;
Zdeltae=coeff*CLiH;
Mdeltae=coeff2*CmiH/i_yy;
XdeltaT=coeff*CDdeltaT;
ZdeltaT=coeff*CLdeltaT;
MdeltaT=coeff2*CMdeltaT/i_yy;

%2.	The matrix B of aircraft

B=[Xdeltae,XdeltaT;
  Zdeltae,ZdeltaT;
  Mdeltae+Zdeltae*M_w_point,MdeltaT+ZdeltaT*M_w_point;
  0,0 ];
disp('B = ')
disp(B)

%eta
eta(1,:)=B(:,1).';
eta(2,:)=B(:,2).';
disp('eta = ')
disp(eta)

%characteristic polynom
coeffPoly = poly(A);       
syms x ;
polynome = vpa(simplify(det(A-x*eye(4))),4);

%eigen values
[k,l]=eig(A);

eigenvect1=l(1,1);
eigenvect2=l(2,2);
eigenvect3=l(3,3);
eigenvect4=l(4,4);

%----------Short-period mode-------------
disp('short-period mode');
zeta_sp=sqrt(1/(1+((imag(eigenvect1)/real(eigenvect1))^2)));
disp(zeta_sp);
omega_sp = -(real(eigenvect1))/zeta_sp;
disp(omega_sp);

%----------Phugoid mode------------
zeta_ph = sqrt(1/(1+((imag(eigenvect3)/real(eigenvect3))^2)));
omega_ph = (real(eigenvect3))/zeta_ph;
disp('rad/s');

%---------Curves of longitudinal motion------
T1=(0:1:500);
T2=(0:0.1:10);

figure (1)

u=k(:,2);
[t,x] = ode45('DATA', [0 15], u);
plot(t,x(:,1),t,x(:,2),t,x(:,3),t,x(:,4))

legend('Axial velocity ','angle of attack','pitch rate ','pitch angle ')
xlabel('Time(s)')
title('Short Period Mode')

figure (2)

u=k(:,3);
[t,x] = ode45('DATA', [0 1000], u);
plot(t,x(:,1),t,x(:,2),t,x(:,3),t,x(:,4)) 


legend('Axial velocity ','angle of attack ','pitch rate ','pitch angle ')
xlabel('Time(s)')
title('Phygoid Mode')

%------------------Transfer function----------------------
syms s ;
format short e

%-----------------Axial velocity---------------------------
num1=vpa([Xdeltae g_0; -Zdeltae/speed s],3);
deno1=vpa([s-X_u g_0; (-Z_u/speed) s],3);
detnum1=det(num1);
detdeno1=det(deno1);

TF_U=vpa(detnum1/detdeno1,3);

tfnum1=[- 0.019671637229293992277234792709351 0.00047189993940525493349228708743358];
tfdeno1=[1 0.056204677798632474150508642196655 -0.040196427499300170986772926174551];

tf_U=tf(tfnum1,tfdeno1)

vpa(ilaplace(TF_U),2)

[z,gain]=zero(tf_U)


%-----------------Pitch angle------------------------------
num4=vpa((Zdeltae/speed)+((X_u*Zdeltae)/speed)-((Z_u*Xdeltae)/speed),3);
deno4=vpa((s^2)-(X_u*s)-((Z_u*g_0)/speed),3);

TF_Theta=num4/deno4;

tfnum4=[-1.26e-4];
tfdeno4=[1 0.0562 0.0402];

tf_Theta=tf(tfnum4,tfdeno4)

vpa(ilaplace(TF_Theta),2)
%[z,gain]=zero(tf_Theta)

% -----------------Angle of attack------------------------------
num2=vpa([Zdeltae/speed -1 ; Mdeltae+((M_w_point*Zdeltae)/speed) s-(M_q+M_w_point)],3) ;
deno2=vpa([s-(Z_w/speed) -1 ; (-(M_w+((M_w_point*Z_w)/speed))) (s-(M_q+M_w_point))],3) ;

detnum2=det(num2);
detdeno2=det(deno2);

tfnum2=[-0.000048103969358059828209661645814776 -0.020865074388508721822862459897003];
tfdeno2=[1 0.3239312194041303882841020822525 0.0082443925494719386795529395837908];

TF_alpha=detnum2/detdeno2;

tf_alpha=tf(tfnum2,tfdeno2)

vpa(ilaplace(TF_alpha),2)

%[z,gain]=zero(tf_alpha)

% -----------------Pitch rate------------------------------
 
num3=vpa([s-(Z_w/speed) Zdeltae/speed ; (-(M_w+((M_w_point*Z_w)/speed))) Mdeltae+((M_w_point*Zdeltae)/speed)],3) ;
deno3=vpa([s-(Z_w/speed) -1 ; (-(M_w+((M_w_point*Z_w)/speed))) (s-(M_q+M_w_point))],3) ;

detnum3=det(num3) ;
detdeno3=det(deno3) ;

tfnum3=[-0.020849781260039890184998512268066 -0.00012506502437473327746899455895901];
tfdeno3=[1 0.3239312194041303882841020822525 0.0082443925494719386795529395837908];
 
TF_q =detnum3/detdeno3;
 
tf_q=tf(tfnum3,tfdeno3)

vpa(ilaplace(TF_q),2)
