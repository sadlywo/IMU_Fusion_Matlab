%% IMU数据融合与姿态解算对比分析
clear; clc; close all;

% 参数设置
fs = 100;            % 采样频率100Hz
dt = 1/fs;           % 采样间隔
t = 0:dt:10;         % 10秒数据
N = length(t);

% 生成真实姿态角（单位：弧度）
% 假设载体做正弦摆动运动 
roll_true = 0.5*sin(2*pi*0.2*t)';    % 转换为列向量
pitch_true = 0.3*sin(2*pi*0.3*t)';   % 转换为列向量
yaw_true = 0.3*rand(length(t),1);        % 列向量,这里的yaw取随机数

% 加速度计真实值（静止状态下测量重力矢量）
accel_true = zeros(N,3);
for i=1:N
    Rx = [1 0 0; 0 cos(roll_true(i)) -sin(roll_true(i)); 0 sin(roll_true(i)) cos(roll_true(i))];
    Ry = [cos(pitch_true(i)) 0 sin(pitch_true(i)); 0 1 0; -sin(pitch_true(i)) 0 cos(pitch_true(i))];
    R = Ry * Rx;
    accel_true(i,:) = (R * [0; 0; 9.81])'; % 重力矢量在载体坐标系下的投影
end

% 陀螺仪真实角速度（通过姿态角微分得到）
gyro_true = zeros(N,3);
gyro_true(:,1) = gradient(roll_true, dt);   % x轴角速度（横滚角变化率）
gyro_true(:,2) = gradient(pitch_true, dt);  % y轴角速度
gyro_true(:,3) = gradient(yaw_true, dt);    % z轴角速度

% 添加传感器误差模型
accel_bias = [0.1 -0.05 0.2];    % 加速度计偏置(m/s²)
gyro_bias = [0.01 -0.02 0.005];  % 陀螺仪偏置(rad/s)
accel_noise = 0.5;              % 加速度计噪声强度 
gyro_noise = 0.05;              % 陀螺仪噪声强度

% 生成带噪声的测量值
accel_meas = accel_true + accel_bias + accel_noise*randn(N,3);
gyro_meas = gyro_true + gyro_bias + gyro_noise*randn(N,3);

%% 原始数据可视化
figure('Position',[100 100 1000 800])

% 加速度计数据对比
subplot(2,1,1)
hold on
plot(t, accel_true(:,1), 'b')
plot(t, accel_meas(:,1), 'b--')
plot(t, accel_true(:,2), 'g')
plot(t, accel_meas(:,2), 'g--')
plot(t, accel_true(:,3), 'r')
plot(t, accel_meas(:,3), 'r--')
title('加速度计数据对比（实线：真实值，虚线：带噪声测量值）')
xlabel('时间 (s)')
ylabel('加速度 (m/s²)')
legend('X轴真实','X轴测量','Y轴真实','Y轴测量','Z轴真实','Z轴测量')

% 陀螺仪数据对比
subplot(2,1,2)
hold on
plot(t, gyro_true(:,1), 'b')
plot(t, gyro_meas(:,1), 'b--')
plot(t, gyro_true(:,2), 'g')
plot(t, gyro_meas(:,2), 'g--')
plot(t, gyro_true(:,3), 'r')
plot(t, gyro_meas(:,3), 'r--')
title('陀螺仪数据对比（实线：真实值，虚线：带噪声测量值）')
xlabel('时间 (s)')
ylabel('角速度 (rad/s)')
legend('X轴真实','X轴测量','Y轴真实','Y轴测量','Z轴真实','Z轴测量')

%% 算法1：基本互补滤波
% 初始化变量
roll_comp = zeros(N,1);
pitch_comp = zeros(N,1);
alpha = 0.98;  % 滤波系数

%% 算法2：Mahony互补滤波
% 初始化变量
roll_mahony = zeros(N,1);
pitch_mahony = zeros(N,1);
ki = 0.01;     % 积分增益
kp = 0.5;      % 比例增益
eInt = [0 0 0];% 误差积分

%% 算法3：Madgwick梯度下降
% 初始化四元数
q = [1 0 0 0]';  % 初始姿态四元数
beta = 0.1;       % 梯度下降步长

%% 算法4：卡尔曼滤波
% 状态变量：[横滚角, 俯仰角, 陀螺x偏差, 陀螺y偏差]
x = [0; 0; 0; 0];
P = eye(4);         % 状态协方差矩阵
Q = diag([0.01, 0.01, 0.001, 0.001]); % 过程噪声
R = diag([0.1, 0.1]);                 % 观测噪声

% 存储结果
roll_kalman = zeros(N,1);
pitch_kalman = zeros(N,1);

for k=2:N
    % 陀螺仪积分
    gyro_rate = [gyro_meas(k,1), gyro_meas(k,2), gyro_meas(k,3)];
    roll_gyro = roll_comp(k-1) + gyro_rate(1)*dt;
    pitch_gyro = pitch_comp(k-1) + gyro_rate(2)*dt;
    
    % 加速度计估计姿态
    ax = accel_meas(k,1);
    ay = accel_meas(k,2); 
    az = accel_meas(k,3);
    roll_acc = atan2(ay, sqrt(ax^2 + az^2));
    pitch_acc = atan2(-ax, sqrt(ay^2 + az^2));
    
    % 互补滤波
    roll_comp(k) = alpha*roll_gyro + (1-alpha)*roll_acc;
    pitch_comp(k) = alpha*pitch_gyro + (1-alpha)*pitch_acc;
end

%% 算法实现
% 基本互补滤波
for k=2:N
    % 原互补滤波实现代码...
end

% Mahony互补滤波
for k=2:N
    % 陀螺仪积分
    gyro_rate = [gyro_meas(k,1), gyro_meas(k,2), gyro_meas(k,3)];
    roll_gyro = roll_mahony(k-1) + gyro_rate(1)*dt;
    pitch_gyro = pitch_mahony(k-1) + gyro_rate(2)*dt;
    
    % 加速度计估计姿态
    ax = accel_meas(k,1); ay = accel_meas(k,2); az = accel_meas(k,3);
    roll_acc = atan2(ay, sqrt(ax^2 + az^2));
    pitch_acc = atan2(-ax, sqrt(ay^2 + az^2));
    
    % 误差修正
    error = [roll_acc - roll_gyro, pitch_acc - pitch_gyro, 0];
    eInt = eInt + error * dt;
    
    % 融合
    roll_mahony(k) = roll_gyro + kp*error(1) + ki*eInt(1);
    pitch_mahony(k) = pitch_gyro + kp*error(2) + ki*eInt(2);
end

% 初始化四元数
q = [1 0 0 0]';  % 初始姿态四元数
beta = 0.1;       % 梯度下降步长

% 存储结果
roll_madg = zeros(N,1);
pitch_madg = zeros(N,1);

% Madgwick算法实现
for k=1:N
    % 陀螺仪测量值
    gyro = gyro_meas(k,:);
    
    % 加速度计归一化
    acc = accel_meas(k,:)/norm(accel_meas(k,:));
    
    % 梯度下降优化
    F = [2*(q(2)*q(4) - q(1)*q(3)) - acc(1);
         2*(q(1)*q(2) + q(3)*q(4)) - acc(2);
         2*(0.5 - q(2)^2 - q(3)^2) - acc(3)];
    J = [-q(3)  q(4)  -q(1)  q(2);
          q(2)  q(1)   q(4)  q(3);
          0    -2*q(2) -2*q(3) 0];
    step = (J'*F);
    step = step/norm(step);
    q_dot = 0.5*qmult(q)*[0 gyro]' - beta*step;
    
    % 四元数更新
    q = q + q_dot*dt;
    q = q/norm(q);
    
    % 转换为欧拉角
    roll_madg(k) = atan2(2*(q(1)*q(2) + q(3)*q(4)), 1 - 2*(q(2)^2 + q(3)^2));
    pitch_madg(k) = asin(2*(q(1)*q(3) - q(4)*q(2)));
end

%% 卡尔曼滤波实现
for k=2:N
    % 预测步骤
    F = [1 0 -dt 0; 
         0 1 0 -dt;
         0 0 1 0;
         0 0 0 1];
    x = F * x;
    P = F * P * F' + Q;
    
    % 更新步骤
    H = [1 0 0 0;
         0 1 0 0];
    z = [atan2(accel_meas(k,2), sqrt(accel_meas(k,1)^2 + accel_meas(k,3)^2));
         atan2(-accel_meas(k,1), sqrt(accel_meas(k,2)^2 + accel_meas(k,3)^2))];
    y = z - H * x;
    S = H * P * H' + R;
    K = P * H' / S;
    
    x = x + K * y;
    P = (eye(4) - K * H) * P;
    
    % 存储结果
    roll_kalman(k) = x(1);
    pitch_kalman(k) = x(2);
end

%% 结果可视化
figure('Position',[100 100 1000 800])

% 横滚角对比
subplot(2,1,1)
plot(t, rad2deg(roll_true), 'k--', 'LineWidth',1.5)
hold on
plot(t, rad2deg(roll_comp), 'b')
plot(t, rad2deg(roll_mahony), 'g')
plot(t, rad2deg(roll_madg), 'r')
plot(t, rad2deg(roll_kalman), 'm')
title('横滚角对比')
xlabel('时间 (s)')
ylabel('角度 (°)')
legend('真实值','基本互补','Mahony','Madgwick','卡尔曼')

% 俯仰角对比
subplot(2,1,2)
plot(t, rad2deg(pitch_true), 'k--', 'LineWidth',1.5)
hold on
plot(t, rad2deg(pitch_comp), 'b')
plot(t, rad2deg(pitch_mahony), 'g')
plot(t, rad2deg(pitch_madg), 'r')
plot(t, rad2deg(pitch_kalman), 'm')
title('俯仰角对比')
xlabel('时间 (s)')
ylabel('角度 (°)')
legend('真实值','基本互补','Mahony','Madgwick','卡尔曼')

%% RMSE计算
rmse_comp_roll = sqrt(mean((roll_comp - roll_true).^2));
rmse_madg_roll = sqrt(mean((roll_madg - roll_true).^2));

rmse_comp_pitch = sqrt(mean((pitch_comp - pitch_true).^2));
rmse_madg_pitch = sqrt(mean((pitch_madg - pitch_true).^2));

% 计算各算法RMSE
rmse_comp_roll = sqrt(mean((roll_comp - roll_true).^2));
rmse_mahony_roll = sqrt(mean((roll_mahony - roll_true).^2));
rmse_madg_roll = sqrt(mean((roll_madg - roll_true).^2));
rmse_kalman_roll = sqrt(mean((roll_kalman - roll_true).^2));

rmse_comp_pitch = sqrt(mean((pitch_comp - pitch_true).^2));
rmse_mahony_pitch = sqrt(mean((pitch_mahony - pitch_true).^2));
rmse_madg_pitch = sqrt(mean((pitch_madg - pitch_true).^2));
rmse_kalman_pitch = sqrt(mean((pitch_kalman - pitch_true).^2));

fprintf('【算法对比】\n');
fprintf('基本互补滤波 - 横滚角: %.4f°, 俯仰角: %.4f°\n',...
    rad2deg(rmse_comp_roll), rad2deg(rmse_comp_pitch));
fprintf('Mahony滤波   - 横滚角: %.4f°, 俯仰角: %.4f°\n',...
    rad2deg(rmse_mahony_roll), rad2deg(rmse_mahony_pitch)); 
fprintf('Madgwick算法 - 横滚角: %.4f°, 俯仰角: %.4f°\n',...
    rad2deg(rmse_madg_roll), rad2deg(rmse_madg_pitch));
fprintf('卡尔曼滤波   - 横滚角: %.4f°, 俯仰角: %.4f°\n',...
    rad2deg(rmse_kalman_roll), rad2deg(rmse_kalman_pitch));

%% 四元数乘法函数
function ab = qmult(a)
    ab = [a(1) -a(2) -a(3) -a(4);
          a(2)  a(1) -a(4)  a(3);
          a(3)  a(4)  a(1) -a(2);
          a(4) -a(3)  a(2)  a(1)];
end
