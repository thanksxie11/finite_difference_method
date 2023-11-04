save('dataFile.mat', 'Psy1x', 'Psx1x', 'net1x', 'Psy2x', 'Psx2x', 'net2x');
[Dec,Obj,Con]=platemo('objFcn',@fP,'conFcn',{@fPs @find_a1},'encoding',[1,1,1,1],...
'lower',[60,2000000,4,7.65],'upper',[60,2000000,6.5,10.5]);
% m=[60 2000000 3.62 10.5;60 2000000 3.62 9.9407;60 2000000 4.62 9.9407;60 2000000 5.62 9.9407;60 2000000 6.62 9.9407]';
% m=[60 2000000 4.3 9.9407];
% p=fP(m);
% p2=fPs(m);
% p3=find_a1(m);

function PopObj = fP(x)
    % global Psy Psx net ;  
    load('dataFile.mat');
    % 对 x 进行归一化
    x_new_norm = mapminmax('apply',[x(3) x(4)]' , Psx1x);
    
    % 使用神经网络进行预测
    y_pred = sim(net1x, x_new_norm);
    
    % 对预测结果进行反归一化
    y = y_pred';
    y_real = mapminmax('reverse', y, Psy1x);
    
    % 取平均值作为目标函数的输出
    PopObj = mean(y_real);
end

function PopObj = fPs(x)
    % global Psy Psx net ;  
    load('dataFile.mat');
    % 对 x 进行归一化
    x_new_norm = mapminmax('apply',[x(1) x(3) x(4)]' , Psx2x);
    
    % 使用神经网络进行预测
    y_pred = sim(net2x, x_new_norm);
    
    % 对预测结果进行反归一化
    y = y_pred';
    y_real = mapminmax('reverse', y, Psy2x);
    
    % 取平均值作为目标函数的输出
    PopObj = mean(y_real)-x(2);
end

function a2 = find_a1(x)
    
    b1 = 0.17; % 请替换为实际值
    a2 = 13; % 请替换为实际值
    b2 = 0.17; % 请替换为实际值
    target_prob = 0.01;
    threshold = 1e-6; % 精度
    
    % 初始的a1搜索范围
    lower_bound = 5;
    upper_bound = 12;

    % 初始化a1为其平均值
    a1 = (lower_bound + upper_bound) / 2;

    % 计算Z的参数 (Z = Y - X)
    mu_Z = a2 - a1;
    sigma_Z = sqrt(b1^2 + b2^2);

    % 计算P(Y - X < 2.5)的概率
    prob = normcdf(2.5, mu_Z, sigma_Z);

    while abs(prob - target_prob) > threshold
        if prob < target_prob
            lower_bound = a1;
        else
            upper_bound = a1;
        end
        a1 = (lower_bound + upper_bound) / 2;

        % 重新计算Z的参数并得到新的概率
        mu_Z = a2 - a1;
        sigma_Z = sqrt(b1^2 + b2^2);
        prob = normcdf(2.5, mu_Z, sigma_Z);
    end
    a2=x(4)-a1;
end

% function P_P=fPs(x)
%     %定义模型参数
%     h2=x(1)/1000000;
%     Ps=500000;
%     b=x(3)/1000;%%封油边宽度
%     R1=0.239;%%油垫封油边内半径，单位m
%     R2=0.239+b;%%油垫油兜内半径，单位m
%     R4=0.308;%%油垫封油边外半径，单位m
%     R3=0.308-b;%%油垫油兜外半径，单位m
%     l=0.01;%%横向封油边宽度
%     phib=l*2*360/(R1+R4)/2/pi;%%油垫角度处的封油边
%     phi2=57.8/360*2*pi;%%油垫扇形角度，前面角度，总体弧度
%     phi1=(57.8-2*phib)/360*2*pi;%%油兜扇形角度，前面角度，总体弧度
%     R_middle=0.5*(R1+R4);%%油垫中心半径
%     phi_middle=0.5*phi2;%%油垫中心弧度
%     pocket_x_scale=(R3-R2)/(R4-R1);%%油兜径向大小
%     pocket_y_scale=phi1/phi2;%%油兜周向大小
%     eta0=0.008;
%     % 
%     A=phi2*(R3*R3-R2*R2)/2;
%     % syms h1 h2 Ps a;
%     lz=0.066;%%中间沟壑长度
%     bz=0.002;%%中间节流缝隙长度
%     Lf=(2*R4-b)*phi2/2-2*b+(2*R1+b)*phi2/2+(2*R4-2*R1-b);%封油边长度
%     ll=0.09;%%节流器两侧缝隙长度，单位m。
%     bl=x(4)/1000;%%节流器两侧缝隙宽度，单位m。
%     lk=0.008;%%节流器两侧缝隙宽处长度，单位m。
%     bk=0.01;%%节流器两侧缝隙宽出宽度，单位m。
%     lm=0.005;%%节流器中间窜油长度，单位m。
%     bm=0.0185;%%节流器中间窜油宽度，单位m。
% 
%     %各个液阻计算值
%     r2c=12*eta0*b/Lf/h2^3;%下油垫出油液阻
%     Rz2=6*eta0*bm/lm/(h2)^3;%下中间油垫窜油液阻
%     Rh22=6*eta0*bz/lz/(h2)^3;%下进油液阻
%     % r2=1/((1/Rh22)+(1/Rz1));%上综合液阻
% 
%     % r1=12*a*bz/h1^3/lz/2;%上油垫进油液阻
%     % r2=12*a*bz/h2^3/lz/2;%下油垫进油液阻
%     Rnl2=6*eta0*bl/ll/(h2)^3;%下窜油液阻
%     Rnk2=3*eta0*bk/lk/(h2)^3;%下窜油宽处液阻
%     r=1/((1/Rnl2)+(1/Rnk2));
% 
%     for ii=1:1:5000
%         Ps=Ps+10000;
%         P01=Ps*(r+r2c)/(r+r2c+Rh22);
%         P02=Ps*r2c/(r+r2c+Rh22);
%         [Wy,Qy]=fWy(P02,Ps,h2,b);
%         [Wj,Qj]=fWj(P02,Ps,P01,h2,b,bl);
%         W02=Wy+Wj;
% 
%         W01=P01*A;
%         % p11(ii)=P01;
%         % p22(ii)=P02;
% 
%         if((W01-W02)>=7200)
%             break;
%         end
% 
%     end
%     %     Ps11(jj)=Ps;
%     %     b=b-0.001;
%     P_P=Ps-x(2);
% end
% 
% function [Wy,Qy]=fWy(Pb,Ps,H,b)
%     b=b;
%     P0=Pb;
%     Ps=Ps;
%     H0=H;
%     Initialising_youdian
%     FDA_youdian
%     Performance_youdian
%     Wy=W;
%     Qy=Q;
% end
% 
% function [Wj,Qj]=fWj(Pb,Ps,Pd,H,b,bl)
%     b=b;
%     P0=Pb;
%     Ps=Ps;
%     Pd=Pd;
%     H0=H;
%     Initialising_jieliuqi
%     FDA_jieliuqi
%     Performance_jieliuqi
%     Wj=W;
%     Qj=Q;
% end
