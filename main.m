clear all
clc

r0=[0.25,0.25];
A = []; % 線性不等式拘束條件的係數矩陣
b = []; % 線性不等式拘束條件的係數向量 AX <= b
Aeq = []; % 線性不等式拘束條件的係數向量
beq = []; % 線性等式拘束條件的係數向量 AeqX = beq
lb = [0.001; 0.001]; % 設計空間的upper bounds
ub = [0.5; 0.5]; % 設計空間的lower bounds
options = optimset ('display','off','Algorithm','sqp');
[r,fval,exitflag] = fmincon(@(r)object_function(r), r0, A, b, Aeq, beq, lb, ub, @(r)nonlcon(r),options);

%繪製設計空間、可行解空間與目標函數值圖
YFT=1;
MDT=1;
YF=zeros(2500,2);
MD=zeros(2500,2);

for (i=1:50)
    for(j=1:50)
       
[sigma,Q]=sol_TenBarTruss(i*0.01,j*0.01);
if (max(sigma)>2.5*10^8 | min(sigma)<-2.5*10^8)
    %figure(1)
   YF(YFT,1)=i*0.01;
   YF(YFT,2)=0.01*j;
   YFT=YFT+1;
    %plot(i*0.01,j*0.01,'gX');
    %hold on
end
if (sqrt(Q(3)^2+Q(4)^2)>0.02) 
   
    %plot(i*0.01,j*0.01,'bX');
    MD(MDT,1)=i*0.01;
    MD(MDT,2)=0.01*j;
    MDT=MDT+1;
    %hold on
end
    end
end

x = 0.01:0.01:0.5;
y = 0.01:0.01:0.5;

[X,Y] = meshgrid(x,y);
Z=6*X.^2*pi()*0.914*7860+4*Y.^2*pi()*7860*0.914*sqrt(2);
figure(1)
contour(X,Y,Z,7);
hold on

plot(YF(1:YFT,1),YF(1:YFT,2),'gX');
hold on
plot(MD(1:MDT,1),MD(1:MDT,2),'bX');
hold on
plot(r(1),r(2),'r*');

axis([0.01,0.5,0.01,0.5]);
xlabel('r1(m)');
ylabel('r2(m)');

legend('Weight','Yielding stress>250*10^6Pa','Max displacement>0.02m','optical point');
