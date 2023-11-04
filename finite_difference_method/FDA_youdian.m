%%%  ���޲���㷨  Finite Difference Algorithm

NUM=500;%%������������
w_sor_p=1.77;%%�ɳ�����ȡֵ��0-2֮��
tolerance=1e-7;%%���㾫��
error=0;%%�������
error_max=0;%%�˴ε��������е�������
p_old=0;%%ǰһ�ε�ѹ������ֵ���ж��Ƿ�������

p=zeros(y_pieces,x_pieces);%%��ʼ��ѹ��Ϊ0
dpdx=zeros(y_pieces,x_pieces-1);%%��ʼ��ѹ���ĵ���Ϊ0
dpdy=zeros(y_pieces,x_pieces-1);%%��ʼ��ѹ���ĵ���Ϊ0
%%v=zeros(y_pieces,x_pieces);%%��ʼ�����ٶ�Ϊ0
v_x=zeros(y_pieces,x_pieces,z_pieces);%%��ʼ��x���ٶ�Ϊ0
v_y=zeros(y_pieces,x_pieces,z_pieces);%%��ʼ��y���ٶ�Ϊ0

%%ѹ���߽�����
for j=1:1:y_pieces %%�Ͷ���ȫ��1����ȫ��0
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
                p(j,i)=((y_step)^2*x(i)*(h(j,i))^3/eta(j,i)*p(j,i+1)+(y_step)^2*x(i-1)*(h(j,i-1))^3/eta(j,i-1)*p(j,i-1)+(x_step)^2*(h(j,i))^3/eta(j,i)/x(i)*p(j+1,i)+(x_step)^2*(h(j-1,i))^3/eta(j-1,i)/x(i)*p(j-1,i)-6*x_step*(y_step)^2*(x(i)*h(j,i)*U_x(j,i)-x(i-1)*h(j,i-1)*U_x(j,i-1))-6*y_step*(x_step)^2*(h(j,i)*U_y(j,i)-h(j-1,i)*U_y(j-1,i))-12*(x_step)^2*(y_step)^2*dhdt-x_step*(y_step)^2*(h(j,i)^3*U_y(j,i)*rho(j,i)/eta(j,i)-h(j,i-1)^3*U_y(j,i-1)*rho(j,i-1)/eta(j,i-1)))/((y_step)^2*x(i)*(h(j,i))^3/eta(j,i)+(y_step)^2*x(i-1)*(h(j,i-1))^3/eta(j,i-1)+(x_step)^2*(h(j,i))^3/x(i)/eta(j,i)+(x_step)^2*(h(j-1,i))^3/x(i)/eta(j-1,i));%%��ֹ�ʽ
                p(j,i)=w_sor_p*p(j,i)+(1-w_sor_p)*p_old;%%��γ��ɳڣ���p��p_old�ļ�Ȩ��
                error=abs(p(j,i)-p_old);%%����һ�ε����
                if error_max<error
                    error_max=error;%%ѡ���˴ε���������������
                end
            end
        end
    end
    drawnow;
    if error_max<tolerance %%���������ж�
       break; 
    end
    if n~=NUM
        error_max=0;%%��һ�ε��������ж�
    end
end

%%����ѹǿ�Ľ����ٶȷֲ�
for j=1:1:y_pieces-1
    for i=1:1:x_pieces-1
        dpdx(j,i)=(p(j,i+1)-p(j,i))/x_step;%%����ѹǿ�ĵ���
        dpdy(j,i)=(p(j+1,i)-p(j,i))/y_step; 
        for k=1:1:z_pieces
            if x(i)<0.5*(x_max+x_min)%%��Ϊ��������һ�����ǰ��ֿ���
                v_x(j,i,k)=((z(k)*h(j,i))^2-z(k)*h(j,i))*0.5*dpdx(j,i)/eta(j,i)+U_x(j,i)*z(k)*h(j,i);
            elseif x(i)>=0.5*(x_max+x_min)%% ��>=�� ����Ϊ�㷨������һλ��ԭ��
                v_x(j,i+1,k)=((z(k)*h(j,i))^2-z(k)*h(j,i))*0.5*dpdx(j,i)/eta(j,i)+U_x(j,i)*z(k)*h(j,i);
            end
            v_x(j,round((x_pieces+1)/2),k)=((z(k)*h(j,round((x_pieces+1)/2)))^2-z(k)*h(j,round((x_pieces+1)/2)))*0.5*dpdx(j,round((x_pieces+1)/2))/eta(j,round((x_pieces+1)/2))+U_x(j,round((x_pieces+1)/2))*z(k)*h(j,round((x_pieces+1)/2));
            if y(j)<0.5*(y_max+y_min)%%��Ϊ��������һ�����ǰ��ֿ���
                v_y(j,i,k)=((z(k)*h(j,i))^2-z(k)*h(j,i))*0.5*dpdy(j,i)/eta(j,i)/x(i)+U_y(j,i)*z(k)*h(j,i);
            elseif y(j)>=0.5*(y_max+y_min)%% ��>=�� ����Ϊ�㷨������һλ��ԭ��
                v_y(j+1,i,k)=((z(k)*h(j,i))^2-z(k)*h(j,i))*0.5*dpdy(j,i)/eta(j,i)/x(i)+U_y(j,i)*z(k)*h(j,i);
            end
            v_y(round((y_pieces+1)/2),i,k)=((z(k)*h(round((y_pieces+1)/2),i))^2-z(k)*h(round((y_pieces+1)/2),i))*0.5*dpdy(round((y_pieces+1)/2),i)/eta(round((y_pieces+1)/2),i)/x(i)+U_y(round((y_pieces+1)/2),i)*z(k)*h(round((y_pieces+1)/2),i);
        end
    end
end
v_y(:,x_pieces,:)=2*v_y(:,x_pieces-1,:)-v_y(:,x_pieces-2,:);%%���һ���������
v_x(y_pieces,:,:)=2*v_x(y_pieces-1,:,:)-v_x(y_pieces-2,:,:);%%���һ���������
v=sqrt((v_x.^2)+(v_y.^2));%%���ٶ�