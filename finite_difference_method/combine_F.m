% clc
% clear
NUM_M=5000;
H000=3e-5;
Ps=2E6;
delth=0;
b=0.0035;
% eta0=0.008;
W_G=7200;%负载3t，轴向6组油垫，平均每组油垫承重为W_G
S=zeros(NUM_M,1);
% omega=0*2*pi/60;%%转速，单位，r/min，转/分
Wz=0;
bl=0.0099407;
% %%定义模型参数
% b=0.01;%%封油边宽度，单位m
% R1=0.239;%%油垫封油边内半径，单位m
% R2=0.239+b;%%油垫油兜内半径，单位m
% R4=0.308;%%油垫封油边外半径，单位m 
% R3=0.308-b;%%油垫油兜外半径，单位m
% % l=0.01;%%横向封油边宽度
% l=b;
% phib=l*2*360/(R1+R4)/2/pi;%%油垫角度处的封油边
% phi2=57.8/360*2*pi;%%油垫扇形角度，前面角度，总体弧度
% phi1=(57.8-2*phib)/360*2*pi;%%油兜扇形角度，前面角度，总体弧度
% R_middle=0.5*(R1+R4);%%油垫中心半径
% phi_middle=0.5*phi2;%%油垫中心弧度
% pocket_x_scale=(R3-R2)/(R4-R1);%%油兜径向大小
% pocket_y_scale=phi1/phi2;%%油兜周向大小

for iii=1:1:NUM_M 
    h1=H000-delth;
    h2=H000+delth;
    u613
%     P01=Ps/(1+((H0-delth)^3/(H0+delth)^3));
%     P02=Ps/(1+((H0+delth)^3/(H0-delth)^3));
    %上油垫承载力
    P0=P01;
    H0=h1;
    Initialising_youdian
    FDA_youdian
    Performance_youdian
    W1=W;
    Pd=P02;
    Initialising_jieliuqi
    FDA_jieliuqi
    Performance_jieliuqi
    W2=W;
    Ws=W1+W2;

    %下油垫承载力
    P0=P02;
    H0=h2;
    Initialising_youdian
    FDA_youdian
    Performance_youdian
    W1=W;
    Pd=P01;
    Initialising_jieliuqi
    FDA_jieliuqi
    Performance_jieliuqi
    W2=W;
    Wx=W1+W2;
    
    %求解刚度   
    S(iii)=(Ws-Wx-Wz)/delth;
    Wz=Ws-Wx;

    if abs(Ws-Wx-W_G)>(0.01*W_G)
        delth=delth+0.01e-6;
    else
         break;

    end

end

%求解刚度


qqq1=P01/r1c;
qqq2=P02/r2c;
Q_20_1=qqq1+qqq2;




% [P01,P02]=fP(eta0, b, h1, h2, Ps);
% 
% 
% 
% function [P01,P02]=fP(eta0, b, h1, h2, Ps)
% %各个液阻计算值
% r1c=12*eta0*b/Lf/h1^3;%上油垫出油液阻
% r2c=12*eta0*b/Lf/h2^3;%下油垫出油液阻
% 
% Rz1=6*eta0*bm/lm/(h1)^3;%上中间油垫窜油液阻
% Rz2=6*eta0*bm/lm/(h2)^3;%下中间油垫窜油液阻
% Rh11=6*eta0*bz/lz/(h1)^3;%上进油液阻
% Rh22=6*eta0*bz/lz/(h2)^3;%下进油液阻
% r1=1/((1/Rh11)+(1/Rz2));%上综合液阻
% r2=1/((1/Rh22)+(1/Rz1));%上综合液阻
% 
% % r1=12*a*bz/h1^3/lz/2;%上油垫进油液阻
% % r2=12*a*bz/h2^3/lz/2;%下油垫进油液阻
% Rnl1=6*eta0*bl/ll/(h1)^3;%上窜油液阻
% Rnl2=6*eta0*bl/ll/(h2)^3;%下窜油液阻
% Rnk1=3*eta0*bk/lk/(h1)^3;%上窜油宽处液阻
% Rnk2=3*eta0*bk/lk/(h2)^3;%下窜油宽处液阻
% Rn=1/((1/Rnl1)+(1/Rnl2)+(1/Rnk1)+(1/Rnk2));
% r=1/((1/Rnl1)+(1/Rnl2)+(1/Rnk1)+(1/Rnk2));
% 
% %液阻简化，流量、压力计算
% ra=r*r2c/(r1c+r2c+r);
% rb=r*r1c/(r1c+r2c+r);
% rc=r1c*r2c/(r1c+r2c+r);
% y=(r1+ra)*(r2+rb)/(r1+ra+r2+rb)+rc;
% 
% Q=Ps/y;
% 
% % Q=q1+q2
% q1=(r1+ra)/(r2+rb)*Q/(1+((r1+ra)/(r2+rb)));
% q2=(r2+rb)/(r1+ra)*Q/(1+((r2+rb)/(r1+ra)));
% P01=Ps-q1*r2;
% P02=Ps-q2*r1;
% end