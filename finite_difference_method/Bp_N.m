clear
clc
% %创建训练样本输入集
% load data;
% train_data =xlsread('data.xlsx', 'A1:D80');  % 输入特征数据
% test_data = readmatrix('sample_data2.csv');
% %归一化
data_x1= xlsread('data2.xlsx', 'A1:A73');
data_x2= xlsread('data2.xlsx', 'C1:D73');
data_xx=[data_x1 data_x2];
data_yy= xlsread('data2.xlsx', 'G1:G73'); 
[data_x,Psx] = mapminmax(data_xx', 0, 1);
[data_y,Psy] = mapminmax(data_yy', 0, 1);
% 拆分数据集
data_x=data_x';
data_y=data_y';


% train_data = data(1:75, :);
% test_data = data(76:80, :);
% data = mapminmax(pred_data', 0, 1)';
% data=data';
% num=2;%对应四个特征
% %建立训练集测试集
% x_train=[data(1:48,1).';data(1:48,2).';data(1:48,3).';data(1:48,4).'];
% x_test=[data(49:51,1).';data(49:51,2).';data(49:51,3).';data(49:51,4).'];
% y_train=data(2:49,num).';
% y_test=data(50:52,num).';
x_train = data_x(1:69, 1:3).';
y_train = data_y(1:69, 1).';
x_test = data_x(70:73, 1:3).';
y_test = data_y(70:73, 1).';
%创建BP神经网络
%创建网络
net=newff(minmax(x_train),[4,1],{'tansig','purelin'},'trainlm');%隐层神经元个数，输出层神经元个数
%设置训练次数
net.trainParam.epochs = 20000;
%设置收敛误差
net.trainParam.goal=0.000001;
%训练网络
[net,tr]=train(net,x_train,y_train);
%在训练集和测试集上的表现
y_train_predict=sim(net,x_train);
y_test_predict=sim(net,x_test);
%作图 分别在训练集和测试集上
subplot(1,2,1)
plot(1:length(y_train_predict),y_train_predict,'*',1:length(y_train_predict),y_train,'o')
title('In Train data')
subplot(1,2,2)
plot(1:4,y_test_predict,'*',1:4,y_test,'o')
title('In Test data')
%求出误差 训练集和测试集
train_error=sum(abs(y_train_predict- y_train))/length(y_train);
test_error=sum(abs(y_test_predict- y_test))/length(y_test);