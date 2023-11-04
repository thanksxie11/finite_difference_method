    % clear
    % clc
    % %创建训练样本输入集
    % load data;
    % train_data =xlsread('data.xlsx', 'A1:D80');  % 输入特征数据
    % test_data = readmatrix('sample_data2.csv');
    % %归一化
    data_xx= xlsread('data5.xlsx', 'C1:D350');
    data_yy= xlsread('data5.xlsx', 'E1:E350'); 
    [data_x,Psx1x] = mapminmax(data_xx', 0, 1);
    [data_y,Psy1x] = mapminmax(data_yy', 0, 1);
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
    x_train = data_x(1:320, 1:2).';
    y_train = data_y(1:320, 1).';
    x_test = data_x(321:350, 1:2).';
    y_test = data_y(321:350, 1).';
    %创建BP神经网络
    %创建网络
    net1x=newff(minmax(x_train),[13,1],{'logsig','purelin'},'trainlm');%隐层神经元个数，输出层神经元个数
    %设置训练次数
    net1x.trainParam.epochs = 20000;
    %设置收敛误差
    net1x.trainParam.goal=0.000001;
    %训练网络
    [net1x,tr]=train(net1x,x_train,y_train);
    %在训练集和测试集上的表现
    y_train_predict=sim(net1x,x_train);
    y_test_predict=sim(net1x,x_test);
    %作图 分别在训练集和测试集上
    subplot(1,2,1)
    plot(1:length(y_train_predict),y_train_predict,'*',1:length(y_train_predict),y_train,'o')
    title('In Train data')
    subplot(1,2,2)
    plot(1:30,y_test_predict,'*',1:30,y_test,'o')
    title('In Test data')
    %求出误差 训练集和测试集
    train_error=sum(abs(y_train_predict- y_train))/length(y_train);
    test_error=sum(abs(y_test_predict- y_test))/length(y_test);