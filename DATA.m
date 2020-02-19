function sys=DATA(t,u)

%,w,q,Theta
Xu= -5.6205e-02; 
Xw= 3.4997e-02;


Zu=-2.8725e-01;
Zw=-4.2154e-01;
Zq=0 ;
Zw_dot =9.7747e-04;

Mu=0;
Mw=-6.3456e-03;
Mq= -3.1578e-01;
Mw_dot= -2.1396e-03;

U0=7.0104e+01;

g=9.81;
theta0=0;


V=u(1);
W=u(2);
Q=u(3);
Theta=u(4);

du =Xu*V+Xw*W +0*Q -g*cos(theta0)*Theta;
dw =Zu/(1-Zw_dot)*V + Zw/(1-Zw_dot)*W+ (Zq+U0)/(1-Zw_dot)*Q+(-g*sin(theta0)/(1-Zw_dot))*Theta;
dq =(Mu+(Mw_dot*Zu)/(1-Zw_dot))*V+ (Mw+(Mw_dot*Zw)/(1-Zw_dot))*W+(Mq+(Mw_dot*(Zq+U0 ))/(1-Zw_dot))*Q-(g*Mw_dot*sin(theta0)/(1-Zw_dot))*Theta;
dtheta = Q ;

sys=[du;dw;dq;dtheta];
end