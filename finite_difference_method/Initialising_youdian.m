%%%  ���������ڼ��������͵��ѹ���ֲ�������˳�� Initialising -> FDA -> Performance or DataShow
%%%  ���̣� d          h^3      dP       d     h^3      dP              d(x*U_x*h)    d(U_y*h)
%%%        ----( x * ------ * ------) + ----(------ * ------)  =  6 * ( --------- + ---------- )   
%%%         dx         eta      dx       dy    eta     x*dy                dx           dy
%%%  x�ǰ뾶��y�ǽǶ�
%%%�����ٻ���
%%%  p=p/P0,x=x/R0,y=y/2/pi,h=h/H0��U_x=U_x/(H^2*P0/R4/eta),U_y=U_y/(H^2*P0/R4/eta),rho=rho/(eta^2*R4^2/H0^4/P0),eta=eta/eta0
%%%  T=T/T0��t=t/(eta*L^2/H^2*P0)


%%����ģ�Ͳ���
% b=0.01;%%���ͱ߿���
% H0=0.000035;%%�����Ĥ��ȣ�m
% P0=580410;%%����ѹ����Pa
% Ps=1500000;


R1=0.239;%%�͵���ͱ��ڰ뾶����λm
R2=0.239+b;%%�͵��Ͷ��ڰ뾶����λm
R4=0.308;%%�͵���ͱ���뾶����λm
R3=0.308-b;%%�͵��Ͷ���뾶����λm
l=b;%%������ͱ߿���
phib=l*2*360/(R1+R4)/2/pi;%%�͵�Ƕȴ��ķ��ͱ�
phi2=57.8/360*2*pi;%%�͵����νǶȣ�ǰ��Ƕȣ����廡��
phi1=(57.8-2*phib)/360*2*pi;%%�Ͷ����νǶȣ�ǰ��Ƕȣ����廡��
R_middle=0.5*(R1+R4);%%�͵����İ뾶
phi_middle=0.5*phi2;%%�͵����Ļ���
pocket_x_scale=(R3-R2)/(R4-R1);%%�Ͷ������С
pocket_y_scale=phi1/phi2;%%�Ͷ������С
%%�������
x_pieces=61;%%x��ָ�ȡ������
y_pieces=81;%%y��ָ�ȡ������
z_pieces=101;%%z��ָ�ȡ������
x_max=1;%%x�᷶Χ�ϱ߽磬�����ٻ�һ��ȡֵΪ1
x_min=R1/R4;%%x�᷶Χ�±߽�
x_step=(x_max-x_min)/(x_pieces-1);%%x�Ჽ�����������������
y_max=phi2;%%y�᷶Χ�ϱ߽�
y_min=0;%%y�᷶Χ�±߽�
y_step=(y_max-y_min)/(y_pieces-1);%%y�Ჽ�����������������
z_max=1;%%y�᷶Χ�ϱ߽磬�����ٻ�һ��ȡֵΪ1
z_min=0;%%y�᷶Χ�±߽�
z_step=(z_max-z_min)/(z_pieces-1);%%y�Ჽ�����������������
%%���幤������
eta0=0.008;%%ճ�ȣ�Pa*s

dhdt=0*1e-5/((H0^3)*P0/eta0/(R4^2));%%�����ǰ�浥λm/s������������
%%������Һ����
rho_0=897.1;%%�ܶȣ�kg/m3
k_cc=0.132;%%����ϵ����w/m/K
c_p=1.926;%%�����ݣ�J/kg/K

x=x_min:x_step:x_max;%%����x��
y=y_min:y_step:y_max;%%����y��
z=z_min:z_step:z_max;%%����z��
h=ones(y_pieces,x_pieces);%%��ʼ����϶����Ϊ1
r=zeros(y_pieces,x_pieces);%%��ʼ�����뺯��Ϊ0,��ת̨����
theta=zeros(y_pieces,x_pieces);%%��ʼ���ǶȺ���Ϊ0,��ת̨��������
eta=ones(y_pieces,x_pieces);%%��ʼ��ճ�Ⱥ���Ϊ1
U_x=zeros(y_pieces,x_pieces);%%��ʼ��x���ٶȺ���
U_y=ones(y_pieces,x_pieces);%%��ʼ��y���ٶȺ���
rho=rho_0/(eta0^2*R4^2/H0^4/P0)*ones(y_pieces,x_pieces);%%��ʼ���ܶȣ�����������

%%����ת��
%%ת��������������һ��
omega=100*2*pi/60;%%ת̨ת�٣���λת/��
for j=1:1:y_pieces
    for i=1:1:x_pieces
        U_y(j,i)=x(i)*R4*omega/(H0^2*P0/R4/eta0);%%������ת��
    end
end

%%������תֱ�����꣬y=0�İ뾶��x��������
xx=zeros(y_pieces,x_pieces);%%ֱ�����������ֵ
yy=zeros(y_pieces,x_pieces);%%ֱ������������ֵ
for j=1:1:y_pieces
    for i=1:1:x_pieces
        xx(j,i)=x(i)*cos(y(j));
        yy(j,i)=x(i)*sin(y(j));
    end
end

delta=0*1e-5/H0/2;%%������б��
phi=2*pi/2;%%��б�׶˶�Ӧ�ĽǶȣ�������
%%�����϶
for j=1:1:y_pieces 
    for i=1:1:x_pieces  
        h(j,i)=1;
    end
end