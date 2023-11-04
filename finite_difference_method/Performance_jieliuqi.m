%%%计算工作时的相关参数，承载力、承载力矩，流量
%%%  无量纲化：w=w/(P0*pi*R4^2)，m_y=m/(P0*pi*R4^3),m_x=m/(P0*pi*R4^3)，q=q/(P0*H^3/eta)

w_vertical=0;%%总承载力，竖直方向
m_center_y=0;%%承载弯矩，对油垫中心取矩,沿y+方向
m_center_x=0;%%承载弯矩，对油垫中心取矩,沿x+方向
q_x=0;%%流量，x方向
q_y=0;%%流量，y方向
friction_x=0;
friction_y=0;
Nf_x=0;
Nf_y=0;
for j=1:1:y_pieces-1
    for i=1:1:x_pieces-1
        w_vertical=w_vertical+0.25*(p(j,i)+p(j+1,i)+p(j,i+1)+p(j+1,i+1))*x_step*y_step*0.5*(x(i+1)+x(i));
        m_center_x=m_center_x+0.25*(p(j,i)+p(j+1,i)+p(j,i+1)+p(j+1,i+1))*x_step*y_step*0.5*(x(i+1)+x(i));
        m_center_y=m_center_y-0.25*(p(j,i)+p(j+1,i)+p(j,i+1)+p(j+1,i+1))*x_step*y_step*0.5*(x(i+1)+x(i));
    end
end

for i=1:1:x_pieces-1
    for k=1:1:z_pieces-1
        q_y=q_y+0.25*(v_y(1,i,k)+v_y(1,i+1,k)+v_y(1,i,k+1)+v_y(1,i+1,k+1))*x_step*0.5*(h(1,i)+h(1,i+1))*z_step;
        q_y=q_y-0.25*(v_y(y_pieces,i,k)+v_y(y_pieces,i+1,k)+v_y(y_pieces,i,k+1)+v_y(y_pieces,i+1,k+1))*x_step*0.5*(h(y_pieces,i)+h(y_pieces,i+1))*z_step;
    end
end
for j=1:1:y_pieces-1
    for k=1:1:z_pieces-1
        q_x=q_x+0.25*(v_x(j,1,k)+v_x(j+1,1,k)+v_x(j,1,k+1)+v_x(j+1,1,k+1))*y_step*x(1)*0.5*(h(j,1)+h(j+1,1))*z_step;
        q_x=q_x-0.25*(v_x(j,x_pieces,k)+v_x(j+1,x_pieces,k)+v_x(j,x_pieces,k+1)+v_x(j+1,x_pieces,k+1))*y_step*x(x_pieces)*0.5*(h(j,x_pieces)+h(j+1,x_pieces))*z_step;
    end
end
% for j=1:1:y_pieces-1
%     for i=1:1:x_pieces-1
%         for k=1:1:z_pieces-1
%             if  x(i)<R2/R4 || x(i)>R3/R4 || y(j)<0.5*(phi2-phi1) || y(j)>0.5*(phi2+phi1) || (x(i)>0.258/R4 && x(i)<0.283/R4 && y(j)>0.5*(phi2-8.93/360*2*pi) && y(j)<0.5*(phi2+8.93/360*2*pi))
%         friction_x=friction_x+eta(j,i)*x_step*y_step*0.5*(x(i+1)+x(i))*(0.25*(v_x(j,i,k)+v_x(j+1,i,k)+v_x(j,i+1,k)+v_x(j+1,i+1,k))-0.25*(v_x(j,i,k+1)+v_x(j+1,i,k+1)+v_x(j,i+1,k+1)+v_x(j+1,i+1,k+1)))/(z_step*h(j,i));
%         friction_y=friction_y+eta(j,i)*x_step*y_step*0.5*(x(i+1)+x(i))*(0.25*(v_y(j,i,k)+v_y(j+1,i,k)+v_y(j,i+1,k)+v_y(j+1,i+1,k))-0.25*(v_y(j,i,k+1)+v_y(j+1,i,k+1)+v_y(j,i+1,k+1)+v_y(j+1,i+1,k+1)))/(z_step*h(j,i));     
%             end
%         end
%     end
% end

% for j=1:1:y_pieces-1
%     for i=1:1:x_pieces-1
%         for k=1:1:z_pieces-1
%             if  x(i)<R2/R4 || x(i)>R3/R4 || y(j)<0.5*(phi2-phi1) || y(j)>0.5*(phi2+phi1) || (x(i)>0.258/R4 && x(i)<0.283/R4 && y(j)>0.5*(phi2-8.93/360*2*pi) && y(j)<0.5*(phi2+8.93/360*2*pi))
%         Nf_x=Nf_x+eta(j,i)*x_step*y_step*0.5*(x(i+1)+x(i))*(0.25*(v_x(j,i,k)+v_x(j+1,i,k)+v_x(j,i+1,k)+v_x(j+1,i+1,k))-0.25*(v_x(j,i,k+1)+v_x(j+1,i,k+1)+v_x(j,i+1,k+1)+v_x(j+1,i+1,k+1)))^2/(z_step*h(j,i));
%         Nf_y=Nf_y+eta(j,i)*x_step*y_step*0.5*(x(i+1)+x(i))*(0.25*(v_y(j,i,k)+v_y(j+1,i,k)+v_y(j,i+1,k)+v_y(j+1,i+1,k))-0.25*(v_y(j,i,k+1)+v_y(j+1,i,k+1)+v_y(j,i+1,k+1)+v_y(j+1,i+1,k+1)))^2/(z_step*h(j,i));     
%             end
%         end
%     end
% end
% Friction=sqrt(friction_x^2+friction_y^2);
% F=Friction*H0*P0*R4;
% Nf=(Nf_x+Nf_y)*H0^3*P0^2/eta0;
q=abs(q_x)+abs(q_y);
W=w_vertical*Ps*R4*yy1;
Q=q*H0^3*P0/eta0;
Nt=2*P0*Q;