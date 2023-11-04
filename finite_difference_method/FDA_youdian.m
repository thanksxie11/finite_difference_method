%%%  有限差分算法  Finite Difference Algorithm

NUM=500;%%迭代次数上限
w_sor_p=1.77;%%松弛因子取值在0-2之间
tolerance=1e-7;%%计算精度
error=0;%%迭代误差
error_max=0;%%此次迭代过程中的最大误差
p_old=0;%%前一次的压力计算值，判断是否收敛用

p=zeros(y_pieces,x_pieces);%%初始化压力为0
dpdx=zeros(y_pieces,x_pieces-1);%%初始化压力的倒数为0
dpdy=zeros(y_pieces,x_pieces-1);%%初始化压力的倒数为0
%%v=zeros(y_pieces,x_pieces);%%初始化总速度为0
v_x=zeros(y_pieces,x_pieces,z_pieces);%%初始化x向速度为0
v_y=zeros(y_pieces,x_pieces,z_pieces);%%初始化y向速度为0

%%压力边界条件
for j=1:1:y_pieces %%油兜内全是1，外全是0
    for i=1:1:x_pieces  
        if x(i)>=R2/R4 && x(i)<=R3/R4  && y(j)>=0.5*(phi2-phi1) && y(j)<=0.5*(phi2+phi1)
            p(j,i)=1;
        end
    end
end


for n=1:1:NUM
    for j=2:1:y_pieces-1
        for i=2:1:x_pieces-1
            if  x(i)<R2/R4 || x(i)>R3/R4 || y(j)<0.5*(phi2-phi1) || y(j)>0.5*(phi2+phi1)
                p_old=p(j,i);
                p(j,i)=((y_step)^2*x(i)*(h(j,i))^3/eta(j,i)*p(j,i+1)+(y_step)^2*x(i-1)*(h(j,i-1))^3/eta(j,i-1)*p(j,i-1)+(x_step)^2*(h(j,i))^3/eta(j,i)/x(i)*p(j+1,i)+(x_step)^2*(h(j-1,i))^3/eta(j-1,i)/x(i)*p(j-1,i)-6*x_step*(y_step)^2*(x(i)*h(j,i)*U_x(j,i)-x(i-1)*h(j,i-1)*U_x(j,i-1))-6*y_step*(x_step)^2*(h(j,i)*U_y(j,i)-h(j-1,i)*U_y(j-1,i))-12*(x_step)^2*(y_step)^2*dhdt-x_step*(y_step)^2*(h(j,i)^3*U_y(j,i)*rho(j,i)/eta(j,i)-h(j,i-1)^3*U_y(j,i-1)*rho(j,i-1)/eta(j,i-1)))/((y_step)^2*x(i)*(h(j,i))^3/eta(j,i)+(y_step)^2*x(i-1)*(h(j,i-1))^3/eta(j,i-1)+(x_step)^2*(h(j,i))^3/x(i)/eta(j,i)+(x_step)^2*(h(j-1,i))^3/x(i)/eta(j-1,i));%%差分公式
                p(j,i)=w_sor_p*p(j,i)+(1-w_sor_p)*p_old;%%逐次超松弛，求p与p_old的加权和
                error=abs(p(j,i)-p_old);%%与上一次的误差
                if error_max<error
                    error_max=error;%%选出此次迭代过程中最大误差
                end
            end
        end
    end
    drawnow;
    if error_max<tolerance %%迭代收敛判断
       break; 
    end
    if n~=NUM
        error_max=0;%%下一次迭代重新判断
    end
end

%%根据压强的解求速度分布
for j=1:1:y_pieces-1
    for i=1:1:x_pieces-1
        dpdx(j,i)=(p(j,i+1)-p(j,i))/x_step;%%先算压强的导数
        dpdy(j,i)=(p(j+1,i)-p(j,i))/y_step; 
        for k=1:1:z_pieces
            if x(i)<0.5*(x_max+x_min)%%因为导数项少一项，所以前后分开算
                v_x(j,i,k)=((z(k)*h(j,i))^2-z(k)*h(j,i))*0.5*dpdx(j,i)/eta(j,i)+U_x(j,i)*z(k)*h(j,i);
            elseif x(i)>=0.5*(x_max+x_min)%% “>=” 是因为算法自身后错一位的原因
                v_x(j,i+1,k)=((z(k)*h(j,i))^2-z(k)*h(j,i))*0.5*dpdx(j,i)/eta(j,i)+U_x(j,i)*z(k)*h(j,i);
            end
            v_x(j,round((x_pieces+1)/2),k)=((z(k)*h(j,round((x_pieces+1)/2)))^2-z(k)*h(j,round((x_pieces+1)/2)))*0.5*dpdx(j,round((x_pieces+1)/2))/eta(j,round((x_pieces+1)/2))+U_x(j,round((x_pieces+1)/2))*z(k)*h(j,round((x_pieces+1)/2));
            if y(j)<0.5*(y_max+y_min)%%因为导数项少一项，所以前后分开算
                v_y(j,i,k)=((z(k)*h(j,i))^2-z(k)*h(j,i))*0.5*dpdy(j,i)/eta(j,i)/x(i)+U_y(j,i)*z(k)*h(j,i);
            elseif y(j)>=0.5*(y_max+y_min)%% “>=” 是因为算法自身后错一位的原因
                v_y(j+1,i,k)=((z(k)*h(j,i))^2-z(k)*h(j,i))*0.5*dpdy(j,i)/eta(j,i)/x(i)+U_y(j,i)*z(k)*h(j,i);
            end
            v_y(round((y_pieces+1)/2),i,k)=((z(k)*h(round((y_pieces+1)/2),i))^2-z(k)*h(round((y_pieces+1)/2),i))*0.5*dpdy(round((y_pieces+1)/2),i)/eta(round((y_pieces+1)/2),i)/x(i)+U_y(round((y_pieces+1)/2),i)*z(k)*h(round((y_pieces+1)/2),i);
        end
    end
end
v_y(:,x_pieces,:)=2*v_y(:,x_pieces-1,:)-v_y(:,x_pieces-2,:);%%最后一项线性拟合
v_x(y_pieces,:,:)=2*v_x(y_pieces-1,:,:)-v_x(y_pieces-2,:,:);%%最后一项线性拟合
v=sqrt((v_x.^2)+(v_y.^2));%%总速度