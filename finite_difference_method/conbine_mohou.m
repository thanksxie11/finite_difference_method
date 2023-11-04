clc;
clear;
%定义模型参数
h1=0.000025;
h2=0.000035;
Ps=1500000;
b=0.01;%%封油边宽度
R1=0.239;%%油垫封油边内半径，单位m
R2=0.239+b;%%油垫油兜内半径，单位m
R4=0.308;%%油垫封油边外半径，单位m
R3=0.308-b;%%油垫油兜外半径，单位m
l=0.01;%%横向封油边宽度
phib=l*2*360/(R1+R4)/2/pi;%%油垫角度处的封油边
phi2=57.8/360*2*pi;%%油垫扇形角度，前面角度，总体弧度
phi1=(57.8-2*phib)/360*2*pi;%%油兜扇形角度，前面角度，总体弧度
R_middle=0.5*(R1+R4);%%油垫中心半径
phi_middle=0.5*phi2;%%油垫中心弧度
pocket_x_scale=(R3-R2)/(R4-R1);%%油兜径向大小
pocket_y_scale=phi1/phi2;%%油兜周向大小
eta0=0.008;
% 

% syms h1 h2 Ps a;
lz=0.066;%%中间沟壑长度
bz=0.002;%%中间节流缝隙长度
Lf=(2*R4-b)*phi2/2-2*b+(2*R1+b)*phi2/2+(2*R4-2*R1-b);%封油边长度
ll=0.09;  %%节流器两侧缝隙长度，单位m。0.9
bl=0.006; %%节流器两侧缝隙宽度，单位m。0.006
lk=0.008; %%节流器两侧缝隙宽处长度，单位m。0.008
bk=0.01;  %%节流器两侧缝隙宽出宽度，单位m。0.01;
lm=0.005; %%节流器中间窜油长度，单位m。0.005
bm=0.0185;%%节流器中间窜油宽度，单位m。0.0185

%各个液阻计算值
r1c=12*eta0*b/Lf/h1^3;%上油垫出油液阻
r2c=12*eta0*b/Lf/h2^3;%下油垫出油液阻

Rz1=6*eta0*bm/lm/(h1)^3;%上中间油垫窜油液阻
Rz2=6*eta0*bm/lm/(h2)^3;%下中间油垫窜油液阻
Rh11=6*eta0*bz/lz/(h1)^3;%上进油液阻
Rh22=6*eta0*bz/lz/(h2)^3;%下进油液阻
r1=1/((1/Rh11)+(1/Rz2));%上综合液阻
r2=1/((1/Rh22)+(1/Rz1));%上综合液阻

% r1=12*a*bz/h1^3/lz/2;%上油垫进油液阻
% r2=12*a*bz/h2^3/lz/2;%下油垫进油液阻
Rnl1=6*eta0*bl/ll/(h1)^3;%上窜油液阻
Rnl2=6*eta0*bl/ll/(h2)^3;%下窜油液阻
Rnk1=3*eta0*bk/lk/(h1)^3;%上窜油宽处液阻
Rnk2=3*eta0*bk/lk/(h2)^3;%下窜油宽处液阻
Rn=1/((1/Rnl1)+(1/Rnl2)+(1/Rnk1)+(1/Rnk2));
r=1/((1/Rnl1)+(1/Rnl2)+(1/Rnk1)+(1/Rnk2));

%液阻简化，流量、压力计算
ra=r*r2c/(r1c+r2c+r);
rb=r*r1c/(r1c+r2c+r);
rc=r1c*r2c/(r1c+r2c+r);
y=(r1+ra)*(r2+rb)/(r1+ra+r2+rb)+rc;

Q=Ps/y;

% Q=q1+q2
q1=(r1+ra)/(r2+rb)*Q/(1+((r1+ra)/(r2+rb)));
q2=(r2+rb)/(r1+ra)*Q/(1+((r2+rb)/(r1+ra)));
P01=Ps-q1*r2;
P02=Ps-q2*r1;
[Wy1,Qy1]=fWy(P01,Ps,h1,b);
[Wj1,Qj1]=fWj(P01,Ps,P02,h1,bl);
W01=Wy1+Wj1;
[Wy2,Qy2]=fWy(P02,Ps,h2,b);
[Wj2,Qj2]=fWj(P02,Ps,P01,h2,bl);
W02=Wy2+Wj2;
W=W01-W02;
function [Wy,Qy]=fWy(Pb,Ps,H,b)
    P0=Pb;
    Ps=Ps;
    H0=H;
    b=b;
    Initialising_youdian
    FDA_youdian
    Performance_youdian
    Wy=W;
    Qy=Q;
end

function [Wj,Qj]=fWj(Pb,Ps,Pd,H,bl)
    P0=Pb;
    Ps=Ps;
    Pd=Pd;
    H0=H;
    bl=bl;
    Initialising_jieliuqi
    FDA_jieliuqi
    Performance_jieliuqi
    Wj=W;
    Qj=Q;
end








