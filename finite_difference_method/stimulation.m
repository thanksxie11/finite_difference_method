clc
clear
NUM_M=5000;
input = xlsread('kkx/data5.xlsx', 'A1:D400'); 
HH=input(:,1)/2000000;
PPs=input(:,2);
delth=0;
bbb=input(:,3)/1000;
blll=input(:,4)/1000;
% eta0=0.008;
W_G=7200;%负载3t，轴向6组油垫，平均每组油垫承重为W_G
E=zeros(401,1);
dt=zeros(401,1);
QQ=zeros(401,1);
for jjj=1:1:400
    b=bbb(jjj);
    H000=HH(jjj);
    Ps=PPs(jjj);
    bl=blll(jjj);
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
    
    
        if ((Ws-Wx)<W_G)
            delth=delth+0.02e-6;
        else
             break;
    
        end
        if (delth>=0.00003)
            break;
        end
    end
    dt(jjj)=delth;
    QQ(jjj)=q1+q2;
    E(jjj)=delth/H000;
    delth=0;

end