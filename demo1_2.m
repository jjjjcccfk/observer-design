%State estimation  and unkonwn input reconstruction via both reduced-order
%and high-order sliding mode observers
clear all;
close all;
clc;
%% Iteration Parameter
iter = 500;
n = 1:iter;
t = n * 0.01;
ts=0.01;
%% System Parameters
A=[0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 1;-1 -5 -10 -10 -5];
B=[0 0;0 0;0 -1;1 0;0 0];
D=[0 0;0 0;0 -1;1 0;0 0];
C=[1 0;0 0;0 0;0 1;0 0]';
Ca=[C(1,:);C(1,:)*A;C(1,:)*A*A;C(2,:)];
d1=2*sin(10*t+3);
d2=-2*sawtooth(2*pi*t);
d=[d1;d2];
u=zeros(2,iter);
x=zeros(5,iter);
y=zeros(2,iter);
%ya1参数
re_v=[zeros(2,1) eye(2); 0 0 0];
ei=[zeros(2,1);1];
Hi=[C(1,:)*B;C(1,:)*A*B;C(1,:)*A*A*B];
ya_1=zeros(3,iter);
ya_1(2,1)=0.2;
ya_1(3,1)=0.3;
% ya1_1hat=zeros(1,iter);
% ya1_2hat=zeros(1,iter);
% ya1_3hat=zeros(1,iter);
% ya1_4hat=zeros(1,iter);
lemma1=16;%蓝线振幅20
lemma2=14.5;%蓝线振荡频率15
lemma3=8;%红线振频5
lemma4=24;%红线频率30
ya1_hat=zeros(4,iter);
ya1_hat(2,1)=0.2;
ya1_hat(3,1)=0.3;
re_v_=[zeros(3,1) eye(3); 0 0 0 0];
Hi_=[C(1,:)*B;C(1,:)*A*B;C(1,:)*A*A*B;0 0];
w_1=zeros(5,iter);
%ya2参数
ya_2=zeros(1,iter);
ya2_hat=zeros(2,iter);
w_2=zeros(3,iter);
%% Iteration of state equation by Euler Equation
for k=1:iter
    x(:,k+1)=x(:,k)+ts*(A*x(:,k)+B*u(:,k)+D*d(:,k));
    y(:,k)=C*x(:,k);
end
%计算y_ai
for k=1:iter
    ya_1(:,k+1)=ya_1(:,k)+ts*(re_v*ya_1(:,k)+ei*C(1,:)*A*A*(A*x(:,k)+D*d(:,k)))+Hi*u(:,k);
end
%计算y_ai_hat
for k=1:iter
   w_1(1,k)=ya1_hat(1,k)-ya_1(1,k);
   w_1(2,k)=lemma1*(abs(w_1(1,k)))^(3/4)*sign(w_1(1,k));
   w_1(3,k)=lemma2*(abs(w_1(2,k)))^(2/3)*sign(w_1(2,k));
   w_1(4,k)=lemma3*(abs(w_1(3,k)))^(1/2)*sign(w_1(3,k));
   w_1(5,k)=lemma4*(abs(w_1(4,k)))^(0/1)*sign(w_1(4,k));
%    ya1_1hat(1,k+1)=ya1_1hat(1,k)+ts*(ya1_2hat(1,k)-w1_1+C(1,:)*B*u(:,k));
%    ya1_1hat(1,k+1)=ya1_2hat(1,k)+ts*(ya1_3hat(1,k)-w1_1+C(1,:)*A*B*u(:,k));
%    ya1_3hat(1,k+1)=ya1_3hat(1,k)+ts*(ya1_4hat(1,k)-w1_1+C(1,:)*A*A*B*u(:,k));
%    ya1_4hat(1,k+1)=-w1_4;
    ya1_hat(:,k+1)=ya1_hat(:,k)+ts*(re_v_*ya1_hat(:,k)+Hi_*u(:,k)-w_1(2:5,k));
end
%计算ya_2
for k=1:iter
   ya_2(:,k+1)=ya_2(:,k)+ts*(C(2,:)*A*x(:,k)+C(2,:)*B*u(:,k)+C(2,:)*D*d(:,k)); 
end
%计算ya2_hat
for k=1:iter
     w_2(1,k)=ya2_hat(1,k)-ya_2(1,k);
     w_2(2,k)=lemma1*(abs(w_2(1,k)))^(1/2)*sign(w_2(1,k));
     w_2(3,k)=lemma2*(abs(w_2(2,k)))^(0/1)*sign(w_2(2,k));
     ya2_hat(1,k+1)=ya2_hat(1,k)+ts*(ya2_hat(1,k)-w_2(2,k)+C(2,:)*B*u(:,k));
     ya2_hat(2,k+1)=ya2_hat(2,k)+ts*(-w_2(3,k));
end
%% 数学分析假定设Qa
% syms Pa
Qa=[45.119 -0.0235 0.006 0.006 0.0353;-0.0235 45.1195 -0.0229 0.0012 0.1730;0.0006 -0.0229 45.1204 -0.0216 0.3458;0.0006 0.0012 -0.0216 45.1215 0.2880;0.0353 0.1730 0.34580 0.2880 80.7585];
La=[0.6887 0.7029 0.6562 2.6215;0.8174 0.8041 1.5049 4.1208;0.0557 2.1738 2.6070 2.0628;6.9433 7.7568 23.8463 93.5256;-1.1731 -5.2344 -9.9902 -9.9343];
% X = lyap(A-La*Ca,Qa);
% norm((A-La*Ca)'*X + X*(A-La*Ca) + Qa);    % 验证解的情况
% Pa=vpasolve((A-La*Ca)'*Pa + Pa*(A-La*Ca)==-Qa);
Fa=[-1.2872 -1.9837 -0.2439 0.3675;-0.2949 -1.0836 -10.60011 0.2439];
Pa=[45.4711 0.2044 0.2949 -1.2872 0.0235;0.2044 44.3448 1.0836 -1.9837 0.0154;0.2949 1.0836 10.6001 -0.2439 0.0000; -1.2872 -1.9837 -0.2439 0.3675 0.0000;0.0235 0.0154 0.0000 0.0000 8.0759];
%% smith正交标准化方程
for i=1:4
    for j=1:4
        A11_(i,j)=A(i,j);
    end
end
for i=1:4
    for j=5:5
        A12_(i,j-4)=A(i,j);
    end
end
for i=5:5
    for j=1:4
        A21_(i-4,j)=A(i,j);
    end
end
for i=5:5
    for j=5:5
        A22_(i-4,j-4)=A(i,j);
    end
end
for i=1:4
    for j=1:4
        P11_(i,j)=Pa(i,j);
    end
end
for i=1:4
    for j=5:5
        P12_(i,j-4)=Pa(i,j);
    end
end
for i=5:5
    for j=1:4
        P21_(i-4,j)=Pa(i,j);
    end
end
for i=5:5
    for j=5:5
        P22_(i-4,j-4)=Pa(i,j);
    end
end
Ka_=inv(P22_)*P21_;
ya_hat=[y(1,1:500);ya1_hat(2:3,1:500);y(2,1:500)];
%% Reduced-order Observer
z2_hat=zeros(1,iter);
x_hat=zeros(5,iter);
for k=1:iter
    z2_hat(1,k+1)=z2_hat(1,k)+ts*((A22_+Ka_*A12_)*z2_hat(1,k)+(Ka_*(A11_-A12_*Ka_)+A21_-A22_*Ka_)*ya_hat(:,k)+[Ka_ 1]*B*u(:,k));
    x_hat(:,k)=[ya_hat(:,k);z2_hat(:,k)-Ka_*ya_hat(:,k)];
end
%% unknown input reconstruction
eplise=[ya1_hat(4,1:500)+C(1,:)*A*A*B*u;ya2_hat(2,1:500)+C(2,:)*B*u];
Ga=[C(1,:)*A*A;C(2,:)];
G=Ga*D;
d_hat=inv(G'*G)*G'*(eplise-Ga*(A*x_hat+B*u));
%% Plot
figure(1)
plot(t,ya1_hat(2,1:500)-ya_1(2,1:500))
hold on
plot(t,ya1_hat(3,1:500)-ya_1(3,1:500),'--')
axis([0 5 -40 40])
legend('the estimated error of ya12','the estimated error of ya13');
title('Fig.1. Estimated error of auxiliary output.')
%plot state error
figure(2)
plot(t,x_hat(1,:)-x(1,1:500));
hold on
plot(t,x_hat(2,:)-x(2,1:500),'--');
plot(t,x_hat(3,:)-x(3,1:500),':.');
plot(t,x_hat(4,:)-x(4,1:500),':');
plot(t,x_hat(5,:)-x(5,1:500),'.');
legend('the estimated error of x1','the estimated error of x2','the estimated error of x3','the estimated error of x4','the estimated error of x5');
axis([0 5 -40 40])
title('Fig.2. State estimated error')
%plot d
figure(3)
subplot(2,1,1)
plot(t,d(1,:))
hold on
plot(t,d_hat(1,:),'--')
axis([0 5 -10 10])
legend('actual d1','estimated d1')
title('Fig.3. Reconstruction of unkonwn input d1')
subplot(2,1,2)
plot(t,d(2,:))
hold on
plot(t,d_hat(2,:),'--')
axis([0 5 -3 3])
legend('actual d2','estimated d2')
title('Fig.4. Reconstruction of unkonwn input d2')



