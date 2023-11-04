clc;
clear;
%定义模型参数
h2=0.00006;
Ps=500000;
b=0.002;%%封油边宽度
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
A=phi2*(R3*R3-R2*R2)/2;
% syms h1 h2 Ps a;
lz=0.066;%%中间沟壑长度
bz=0.002;%%中间节流缝隙长度
Lf=(2*R4-b)*phi2/2-2*b+(2*R1+b)*phi2/2+(2*R4-2*R1-b);%封油边长度
ll=0.09;%%节流器两侧缝隙长度，单位m。
bl=0.006;%%节流器两侧缝隙宽度，单位m。
lk=0.008;%%节流器两侧缝隙宽处长度，单位m。
bk=0.01;%%节流器两侧缝隙宽出宽度，单位m。
lm=0.005;%%节流器中间窜油长度，单位m。
bm=0.0185;%%节流器中间窜油宽度，单位m。

%各个液阻计算值
r2c=12*eta0*b/Lf/h2^3;%下油垫出油液阻
Rz2=6*eta0*bm/lm/(h2)^3;%下中间油垫窜油液阻
Rh22=6*eta0*bz/lz/(h2)^3;%下进油液阻
% r2=1/((1/Rh22)+(1/Rz1));%上综合液阻

% r1=12*a*bz/h1^3/lz/2;%上油垫进油液阻
% r2=12*a*bz/h2^3/lz/2;%下油垫进油液阻
Rnl2=6*eta0*bl/ll/(h2)^3;%下窜油液阻
Rnk2=3*eta0*bk/lk/(h2)^3;%下窜油宽处液阻
r=1/((1/Rnl2)+(1/Rnk2));

%液阻简化，流量、压力计算
% ra=r*r2c/(r1c+r2c+r);
% rb=r*r1c/(r1c+r2c+r);
% rc=r1c*r2c/(r1c+r2c+r);
% y=(r1+ra)*(r2+rb)/(r1+ra+r2+rb)+rc;

% Q=Ps/(r+r2c+Rh22);
% 
% % Q=q1+q2
% q1=(r1+ra)/(r2+rb)*Q/(1+((r1+ra)/(r2+rb)));
% q2=(r2+rb)/(r1+ra)*Q/(1+((r2+rb)/(r1+ra)));
% Ps11=zeros(15);
% % p11=zeros(5000);
% % p22=zeros(5000);
% for jj=1:1:7
    for ii=1:1:5000
        Ps=Ps+10000;
        P01=Ps*(r+r2c)/(r+r2c+Rh22);
        P02=Ps*r2c/(r+r2c+Rh22);
        [Wy,Qy]=fWy(P02,Ps,h2,b);
        [Wj,Qj]=fWj(P02,Ps,P01,h2,b,bl);
        W02=Wy+Wj;
    
        W01=P01*A;
        % p11(ii)=P01;
        % p22(ii)=P02;
    
        if((W01-W02)>=7200)
            break;
        end
    
    end
%     Ps11(jj)=Ps;
%     b=b-0.001;
%     Ps=500000;
% end
function [Wy,Qy]=fWy(Pb,Ps,H,b)
    b=b;
    P0=Pb;
    Ps=Ps;
    H0=H;
    Initialising_youdian
    FDA_youdian
    Performance_youdian
    Wy=W;
    Qy=Q;
end

function [Wj,Qj]=fWj(Pb,Ps,Pd,H,b,bl)
    b=b;
    P0=Pb;
    Ps=Ps;
    Pd=Pd;
    H0=H;
    Initialising_jieliuqi
    FDA_jieliuqi
    Performance_jieliuqi
    Wj=W;
    Qj=Q;
end








