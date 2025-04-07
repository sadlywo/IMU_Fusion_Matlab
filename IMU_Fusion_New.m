%% IMU数据融合与姿态估计
clear; clc; close all;

%% 1. 读取IMU数据
data = load('OUTBACK.txt');
acc = data(:,1:3);    % 加速度计数据 (m/s²)
gyro = data(:,4:6);   % 陀螺仪数据 (rad/s)
dt = 1/60;            % 采样间隔
n = size(data,1);     % 数据点数
time = (0:n-1)*dt;    % 时间向量

%% 2. 初始化变量
% 互补滤波参数
alpha = 0.98;  

% Madgwick滤波参数
beta = 0.1;     

% 初始化四元数
q = [1 0 0 0];  

% 存储结果
euler_comp = zeros(n,3);  % 互补滤波欧拉角
euler_madg = zeros(n,3);  % Madgwick滤波欧拉角
euler_gyro = zeros(n,3);  % 陀螺积分欧拉角

%% 3. 主循环处理
for i = 1:n
    %% 3.1 互补滤波
    % 计算加速度计姿态角 (俯仰和横滚)
    pitch_acc = atan2(acc(i,2), sqrt(acc(i,1)^2 + acc(i,3)^2));
    roll_acc = atan2(-acc(i,1), acc(i,3));
    
    % 陀螺仪积分
    if i == 1
        pitch_gyro = pitch_acc;
        roll_gyro = roll_acc;
        yaw_gyro = 0;
    else
        pitch_gyro = euler_comp(i-1,1) + gyro(i,2)*dt;
        roll_gyro = euler_comp(i-1,2) + gyro(i,1)*dt;
        yaw_gyro = euler_comp(i-1,3) + gyro(i,3)*dt;
    end
    
    % 互补滤波融合
    euler_comp(i,1) = alpha*pitch_gyro + (1-alpha)*pitch_acc;
    euler_comp(i,2) = alpha*roll_gyro + (1-alpha)*roll_acc;
    euler_comp(i,3) = yaw_gyro; % 偏航角仅用陀螺仪积分

    %% 3.2 Madgwick滤波
    q = madgwickUpdate(q, acc(i,:), gyro(i,:), dt, beta);
    euler_madg(i,:) = quat2eul(q, 'ZYX');

    %% 3.3 陀螺仪积分
    if i == 1
        euler_gyro(i,1) = pitch_acc;
        euler_gyro(i,2) = roll_acc;
        euler_gyro(i,3) = 0;
    else
        euler_gyro(i,1) = euler_gyro(i-1,1) + gyro(i,2)*dt;
        euler_gyro(i,2) = euler_gyro(i-1,2) + gyro(i,1)*dt;
        euler_gyro(i,3) = euler_gyro(i-1,3) + gyro(i,3)*dt;
    end
end

%% 4. 转换为角度制
euler_comp = rad2deg(euler_comp);
euler_madg = rad2deg(euler_madg);
euler_gyro = rad2deg(euler_gyro);

%% 5. 可视化结果
figure('Name','姿态角对比','NumberTitle','off')
titles = {'俯仰角(Pitch)', '横滚角(Roll)', '偏航角(Yaw)'};

for k = 1:3
    subplot(3,1,k)
    plot(time, euler_gyro(:,k), 'k--', 'LineWidth',1.5)
    hold on
    plot(time, euler_comp(:,k), 'b', 'LineWidth',1)
    plot(time, euler_madg(:,k), 'r', 'LineWidth',1)
    ylabel('角度 (°)')
    title(titles{k})
    legend('陀螺积分','互补滤波','Madgwick滤波')
    grid on
end
xlabel('时间 (s)')

%% 6. 计算并显示RMSE
rmse_comp = sqrt(mean((euler_gyro - euler_comp).^2));
rmse_madg = sqrt(mean((euler_gyro - euler_madg).^2));

fprintf('互补滤波RMSE误差: Pitch=%.2f°, Roll=%.2f°, Yaw=%.2f°\n', rmse_comp)
fprintf('Madgwick滤波RMSE误差: Pitch=%.2f°, Roll=%.2f°, Yaw=%.2f°\n', rmse_madg)

%% Madgwick滤波更新函数
function q = madgwickUpdate(q, acc, gyro, dt, beta)
    % 标准化加速度计测量值
    acc = acc / norm(acc);
    
    % 梯度下降算法优化
    F = [2*(q(2)*q(4) - q(1)*q(3)) - acc(1)
         2*(q(1)*q(2) + q(3)*q(4)) - acc(2)
         2*(0.5 - q(2)^2 - q(3)^2) - acc(3)];
    
    J = [-2*q(3),  2*q(4), -2*q(1), 2*q(2)
          2*q(2),  2*q(1),  2*q(4), 2*q(3)
               0, -4*q(2), -4*q(3),      0];
    
    step = (J'*F);
    step = step / norm(step);
    
    % 手动实现四元数乘法 (替换quatmultiply)
    w = 0; x = gyro(1); y = gyro(2); z = gyro(3);
    q_gyro = [w, x, y, z];
    q_mult = [q(1)*q_gyro(1) - q(2)*q_gyro(2) - q(3)*q_gyro(3) - q(4)*q_gyro(4),
              q(1)*q_gyro(2) + q(2)*q_gyro(1) + q(3)*q_gyro(4) - q(4)*q_gyro(3),
              q(1)*q_gyro(3) - q(2)*q_gyro(4) + q(3)*q_gyro(1) + q(4)*q_gyro(2),
              q(1)*q_gyro(4) + q(2)*q_gyro(3) - q(3)*q_gyro(2) + q(4)*q_gyro(1)];
    
    % 计算四元数导数
    qDot = 0.5 * q_mult - beta * step';
    
    % 积分四元数
    q = q + qDot * dt;
    q = q / norm(q);
end

%% 四元数转欧拉角函数
function euler = quat2eul(q, sequence)
    % 手动实现四元数转欧拉角(ZYX顺序)
    q = q / norm(q); % 替换quatnormalize
    
    % 提取四元数分量
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);
    
    % 计算欧拉角
    roll = atan2(2*(q0*q1 + q2*q3), 1 - 2*(q1^2 + q2^2));
    pitch = asin(2*(q0*q2 - q3*q1));
    yaw = atan2(2*(q0*q3 + q1*q2), 1 - 2*(q2^2 + q3^2));
    
    euler = [yaw, pitch, roll]; % ZYX顺序
end
