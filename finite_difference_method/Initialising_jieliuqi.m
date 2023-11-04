%%%  本程序用于计算扇形油垫的压力分布，运行顺序： Initialising -> FDA -> Performance or DataShow
%%%  方程： d          h^3      dP       d     h^3      dP              d(x*U_x*h)    d(U_y*h)
%%%        ----( x * ------ * ------) + ----(------ * ------)  =  6 * ( --------- + ---------- )   
%%%         dx         eta      dx       dy    eta     x*dy                dx           dy
%%%  x是半径，y是角度
%%%无量纲化：
%%%  p=p/P0,x=x/R0,y=y/2/pi,h=h/H0，U_x=U_x/(H^2*P0/R4/eta),U_y=U_y/(H^2*P0/R4/eta),rho=rho/(eta^2*R4^2/H0^4/P0),eta=eta/eta0
%%%  T=T/T0，t=t/(eta*L^2/H^2*P0)

% H0=0.000035;%%设计油膜厚度，m
% Pd=926560;%%本侧压力，Pa
% P0=580410;%%对侧压力，Pa
% Ps=1500000;%%内反馈供油压力
% bl=0.006;
%%定义模型参数
R4=0.025;%%b油垫宽度，单位m
b1=bl;
bb=(R4-0.003)/2;%%封油边宽度
R1=0;%%油垫封油边内半径，单位m
R2=R1+bb;%%油垫油兜内半径，单位m
R3=R4-bb;%%油垫油兜外半径，单位m
yy1=0.1;
yy2=0.066;
l=(yy1-yy2)/2;%%横向封油边宽度
pocket_x_scale=(R3-R2)/(R4-R1);%%油兜x大小
pocket_y_scale=yy2/yy1;%%油兜y大小
%%定义参数
x_pieces=26;%%x轴分割取点数量
y_pieces=51;%%y轴分割取点数量
z_pieces=101;%%z轴分割取点数量
x_max=1;%%x轴范围上边界，无量纲化一般取值为1
x_min=0;%%x轴范围下边界
x_step=(x_max-x_min)/(x_pieces-1);%%x轴步长，根据条件计算得
y_max=1;%%y轴范围上边界
y_min=0;%%y轴范围下边界
y_step=(y_max-y_min)/(y_pieces-1);%%y轴步长，根据条件计算得
z_max=1;%%y轴范围上边界，无量纲化一般取值为1
z_min=0;%%y轴范围下边界
z_step=(z_max-z_min)/(z_pieces-1);%%y轴步长，根据条件计算得
%%定义工况参数
eta0=0.008;%%粘度，Pa*s

dhdt=0*1e-5/((H0^3)*P0/eta0/(R4^2));%%阻尼项，前面单位m/s，整体无量纲
%%定义油液参数
rho_0=897.1;%%密度，kg/m3
k_cc=0.132;%%传热系数，w/m/K
c_p=1.926;%%比热容，J/kg/K

x=x_min:x_step:x_max;%%定义x轴
y=y_min:y_step:y_max;%%定义y轴
z=z_min:z_step:z_max;%%定义z轴
h=ones(y_pieces,x_pieces);%%初始化间隙函数为1
r=zeros(y_pieces,x_pieces);%%初始化距离函数为0,到转台中心
theta=zeros(y_pieces,x_pieces);%%初始化角度函数为0,到转台中心连线
eta=ones(y_pieces,x_pieces);%%初始化粘度函数为1
U_x=zeros(y_pieces,x_pieces);%%初始化x向速度函数
U_y=zeros(y_pieces,x_pieces);%%初始化y向速度函数
rho=rho_0/(eta0^2*R4^2/H0^4/P0)*ones(y_pieces,x_pieces);%%初始化密度，算离心力用

%%定义转速
%%转速正负计算结果不一致
omega=0*2*pi/60;%%转台转速，单位转/分
for j=1:1:y_pieces
    for i=1:1:x_pieces
        U_y(j,i)=x(i)*R4*omega/(H0^2*P0/R4/eta0);%%无量纲转速
    end
end

% %%极坐标转直角坐标，y=0的半径在x轴正方向
% xx=zeros(y_pieces,x_pieces);%%直角坐标横坐标值
% yy=zeros(y_pieces,x_pieces);%%直角坐标纵坐标值
% for j=1:1:y_pieces
%     for i=1:1:x_pieces
%         xx(j,i)=x(i)*cos(y(j));
%         yy(j,i)=x(i)*sin(y(j));
%     end
% end

delta=0*1e-5/H0/2;%%无量纲斜率
phi=2*pi/2;%%倾斜底端对应的角度，弧度制
%%定义间隙
for j=1:1:y_pieces 
    for i=1:1:x_pieces  
        h(j,i)=1;
    end
end
