%unkown input reconstruction via interval obsever and state and unkown
%input compension feedback controller designs
clear all;
close all;
clc;
%% Iteration Parameter
iter=1500;
n=1:iter;
t=n*0.01;
%% 原始状态方程
A=[0.25 0 0 0 0.25;0.05 0.35 0.05 0.05 0.05;0.1 0.2 0.35 0.1 -0.15;0.05 0.1 0.05 0.3 -0.2;0 0 0 0 0.5];%状态state n=5
B=[3 1;1 0;0 0;0 0;0 0];%u的阶数m=2
D=[3 1;1 0;0 0;0 0;0 0];%d1，d2，阶数q=2
C=[1 1 1 0 0;0 1 1 0 0;0 0 1 0 0];%y的阶数P=3
F=[1;0;1];%w=0.5sin(k) r=1;
%input and unkown input measurement noise
u1=sin(t);
u2=cos(t);
d1=2*sin(0.8*t+2);
d2=0.5*cos(t+1);
w=0.5*sin(t);
u=[u1;u2];
d=[d1;d2];
%actual state and actual output
x=zeros(5,iter);
y=zeros(3,iter);
%% 扩展矩阵
Ee=[eye(5) zeros(5,1)];%扩展矩阵cell Ee
Ce=[C F];
Ae=[A zeros(5,1)];
R=[Ee;Ce];
M=inv(R'*R);
G=M*[eye(5);zeros(1,5)];
re_v=M*[C';F'];
%% Smith正交化，标准状态方程
T=[0.5 3.5 3 0 0 0;0 -11 -11 0 0 0;0 0.44721 0.89443 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0.44721 0 0 0 0 0.89443];
S=[0.5 3 -0.5;0 -11 0;0.70711 0 0.70711];
A_=T*(G*Ae)*inv(T);
B_=T*G*B;
D_=T*G*D;
C_=S*Ce*inv(T);
for i=1:2
    for j=1:2
        A11_(i,j)=A_(i,j);
    end
end
for i=3:6
    for j=1:2
        A21_((i-2),j)=A_(i,j);
    end
end
for i=1:2
    for j=3:6
        A12_(i,(j-2))=A_(i,j);
    end
end
for i=3:6
    for j=3:6
        A22_((i-2),(j-2))=A_(i,j);
    end
end
for i=1:2
    for j=1:2
        B1_(i,j)=B_(i,j);
    end
end
for i=3:6
    for j=1:2
        B2_((i-2),j)=B_(i,j);
    end
end
for i=3:3
    for j=3:6
        C22_(i-2,j-2)=C_(i,j);
    end
end
pi1_=[eye(2) zeros(2,4)]*T*G*Ae*re_v;
pi2_=[zeros(4,2) eye(4)]*T*G*Ae*re_v;
%标准矩阵的state及state_hat Initialization
z1_=zeros(2,iter);
z2_=zeros(4,iter);
z1_hat=zeros(2,iter);
z2_hat=zeros(4,iter);
y1_=zeros(2,iter);
y2_=zeros(1,iter);
xe_hat=zeros(6,iter);
%% State Equation
%原始方程迭代
for k=1:iter
    x(:,k+1)=A*x(:,k)+B*u(:,k)+D*d(:,k);
    y(:,k)=C*x(:,k)+F*w(k);
end
%% Iterative Observation
%标准方程迭代
for k=1:iter
     z1_(:,k+1)=A11_*z1_(:,k)+A12_*z2_(:,k)+B1_*u(:,k)+pi1_*y(:,k)+d(:,k);
     z2_(:,k+1)=A22_*z2_(:,k)+A21_*z1_(:,k)+B2_*u(:,k)+pi2_*y(:,k);
     y1_(:,k)=z1_(:,k);
     y2_(:,k)=C22_*z2_(:,k);
end
y1_(:,iter+1)=z1_(:,iter+1);
%降维obsever
%判断A22_+LC22_为Schur Matrix
L22_=[0;0;0;0];
eig(A22_+L22_*C22_)
for k=1:iter
    z2_hat(:,k+1)=(A22_+L22_*C22_)*z2_hat(:,k)+A21_*y1_(:,k)+B2_*u(:,k)+pi2_*y(:,k)-L22_*y2_(:,k);
    xe_hat(:,k+1)=inv(T)*[y1_(:,k);z2_hat(:,k)]+re_v*y(:,k);
end
%% Interval observer for y1_
y1_h=zeros(2,1500);
y1_l=zeros(2,1500);
d_h=[max(d1);max(d2)];
d_l=[min(d1);min(d2)];
I_=[0.25 0;0 0.2];
Ls=A11_-I_;
Qe=[eye(2) zeros(2,4)]*T*(eye(6)-re_v*Ce);
xe_m=[max(x(1,:));max(x(2,:));max(x(3,:));max(x(4,:));max(x(5,:));max(w)];
xe_l=[min(x(1,:));min(x(2,:));min(x(3,:));min(x(4,:));min(x(5,:));min(w)];
Qe_m=max(Qe,0);
Qe_l=min(0,-Qe);
y1_h(:,1)=Qe_m*xe_m-Qe_l*xe_l;
y1_l(:,1)=Qe_m*xe_l-Qe_l*xe_m;
for k=1:iter
   y1_h(:,k+1)=A11_*y1_h(:,k)+A12_*z2_hat(:,k)+B1_*u(:,k)+pi1_*y(:,k)+Ls*(y1_(:,k)-y1_h(:,k))+d_h;
   y1_l(:,k+1)=A11_*y1_l(:,k)+A12_*z2_hat(:,k)+B1_*u(:,k)+pi1_*y(:,k)+Ls*(y1_(:,k)-y1_l(:,k))+d_l;
   omeg=diag(y1_h(:,k)-y1_l(:,k));
end
% figure(1)
% plot(t,y1_(1,:));
% hold on
% plot(t,y1_h(1,1:1500));
% plot(t,y1_l(1,1:1500));
%% Unknown Input Reconstruction
alfa=zeros(2,iter+1);
for k=1:(iter+1)
    alfa(:,k)=(pinv(diag(y1_h(:,k)-y1_l(:,k)))*(y1_(:,k)-y1_l(:,k))); 
end
d_hat=zeros(2,iter);
for k=1:iter
   d_hat(:,k)=diag(alfa(:,k+1))*(A11_-Ls)*y1_h(:,k)+(eye(2)-diag(alfa(:,k+1)))*(A11_-Ls)*y1_l(:,k)+(Ls-A11_)*y1_(:,k)+(eye(2)-diag(alfa(:,k+1)))*d_l+diag(alfa(:,k+1))*d_h;
end
%% Plot
%actual state and estimated state
figure(1)
subplot(3,2,1)
plot(t,x(1,1:1500));
hold on
plot(t,xe_hat(1,1:1500),'--');
legend('\fontname{Time New Roman}actual x1(k)','estimated x1(k)');
axis([0 15 -15 15]);
title('\fontname{Time New Roman}The estimations of the states x1')
subplot(3,2,2)
plot(t,x(2,1:1500));
hold on
plot(t,xe_hat(2,1:1500),'--');
legend('\fontname{Time New Roman}actual x2(k)','estimated x2(k)');
axis([0 15 -15 15]);
title('\fontname{Time New Roman}The estimations of the states x2')
subplot(3,2,3)
plot(t,x(3,1:1500));
hold on
plot(t,xe_hat(3,1:1500),'--');
legend('\fontname{Time New Roman}actual x3(k)','estimated x3(k)');
axis([0 15 -5 5]);
title('\fontname{Time New Roman}The estimations of the states x3')
subplot(3,2,4)
plot(t,x(4,1:1500));
hold on
plot(t,xe_hat(4,1:1500),'--');
legend('\fontname{Time New Roman}actual x4(k)','estimated x4(k)');
axis([0 15 -5 5]);
title('\fontname{Time New Roman}The estimations of the states x4')
subplot(3,2,5)
plot(t,x(5,1:1500));
hold on
plot(t,xe_hat(5,1:1500),'--');
legend('\fontname{Time New Roman}actual x5(k)','estimated x5(k)');
axis([0 15 -2 2]);
title('\fontname{Time New Roman}The estimations of the states x5')
subplot(3,2,6)
plot(t,w(1,1:1500));
hold on
plot(t,xe_hat(6,1:1500),'--');
legend('\fontname{Time New Roman}actual w(k)','estimated w(k)');
axis([0 15 -2 2]);
title('\fontname{Time New Roman}The estimations of measurement noise')
%interval observer for y1_
figure(2)
subplot(2,1,1)
plot(t,y1_h(1,1:1500),'--')
hold on
plot(t,y1_(1,1:1500))
plot(t,y1_l(1,1:1500),':')
legend('upper y1,1','actual y1,1','lower y1,1')
axis([0 15 -30 30])
title('\fontname{Time New Roman}(a)')
subplot(2,1,2)
plot(t,y1_h(2,1:1500),'--')
hold on
plot(t,y1_(2,1:1500))
plot(t,y1_l(2,1:1500),':')
legend('upper y1,2','actual y1,2','lower y1,2')
axis([0 15 -30 30])
title('\fontname{Time New Roman}(b)')
%signals of alfa_1 and alfa_2 and unkonwn input reconstruction
figure(3)
subplot(2,2,1)
plot(t,alfa(1,1:1500))
axis([0 15 0 1.5])
title('\fontname{Time New Roman}the signals of {\alpha}_1')
subplot(2,2,2)
plot(t,alfa(2,1:1500))
axis([0 15 0 1.5])
title('\fontname{Time New Roman}the signals of {\alpha}_2')
subplot(2,2,3)
plot(t,d(1,1:1500))
hold on
plot(t,d_hat(1,1:1500),'--')
axis([0 15 -3 3])
legend('actual d1(k)','estimated d1(k)')
title('the unknown input reconstruction of d1(k)')
subplot(2,2,4)
plot(t,d(2,1:1500))
hold on
plot(t,d_hat(2,1:1500),'--')
axis([0 15 -3 3])
legend('actual d2(k)','estimated d2(k)')
title('the unknown input reconstruction of d2(k)')
