%
% 千岛湖湖试数据处理3
% 声线修正：遍历角度，得到时延误差范围内的水平传播距离,等梯度声线跟踪
% 定位算法：高斯牛顿迭代法
%
clc,clear,close all

load sound_speed.mat
load boat_position.mat
load data_position.mat
load time_delay.mat


%% 预处理
tic
% 实验数据导入
z_target1 = abs(data_position(3)); % 深度绝对值
z_target2 = data_position(3); % 深度
deepth_a = sound_speed(1,:); % 声速剖面：深度
soundspeed = sound_speed(2,:); % 声速剖面：声速

time_delay = time_delay/2; % 时延信息

pos_num = size(boat_position,1); % 观测点数量

% 参考点数据
real_distance = sqrt((boat_position(:,1) - data_position(1)).^2 + ...
                     (boat_position(:,2) - data_position(2)).^2 + ...
                     (boat_position(:,3)-data_position(3)).^2);
% 最后一层索引
size_end = find(deepth_a>z_target1,1);% 截止层数
soundspeed_targ = sound_speed(2,1:size_end);% 各层声速
delta_z = deepth_a(2:size_end) - deepth_a(1:size_end - 1); % 层厚
g = (soundspeed(2:size_end) - soundspeed(1:size_end-1))./delta_z;% 各层声速梯度


%% 最小二乘法确定基准点
distance1 = mean(soundspeed_targ).*time_delay;
% 选择第一个基站作为参考
x1 = boat_position(1, 1);
y1 = boat_position(1, 2);
z1 = boat_position(1, 3);
d1_sq = distance1(1); % 注意这里使用平方距离

% 构建线性方程组 Ax = b
A = [];
b = [];
valid_distances = []; % 记录用于权重计算的测量距离

for i = 2:pos_num
    xi = boat_position(i, 1);
    yi = boat_position(i, 2);
    zi = boat_position(i, 3);
    di_sq = distance1(i);      % 第i个距离平方值
    
    % 计算方程系数
    ai = -2*(xi - x1);
    bi = -2*(yi - y1);
    ci = -2*(zi - z1);
    di = di_sq - d1_sq - (xi^2 + yi^2 + zi^2 - x1^2 - y1^2 - z1^2);
    
    A = [A; ai, bi, ci];
    b = [b; di];
    valid_distances = [valid_distances; sqrt(di_sq)]; % 存储实际距离值
end

zero_pos = (A' * A) \ (A' * b); % 核心计算公式
zero_pos(3) = 0;

boat_position(:,1) = boat_position(:,1) - zero_pos(1);
boat_position(:,2) = boat_position(:,2) - zero_pos(2);
data_position = data_position - zero_pos';


%% 分层常梯度声线跟踪法（ 遍历入射角确定水平位移 ）
% 初始化参数
sum_y = [];
sum_t = [];
p_exm = [];

% 遍历角度计算传播时间与水平距离 此处入射角为与法平面的夹角
delangal = pi/1800; % 角度分辨率
angal = pi/180:delangal:pi/2;
angal_i = angal;
p = sin(angal)./soundspeed_targ(1); % 折射定律长度
delta_tmin = 1e-4; % 时间精度


for i = 1 : pos_num
    
    % 第一层处理
    size_start = find(deepth_a>boat_position(i,3),1)-1;% 声波发出层数
    g_i = g(size_start:size_end-1);% 各层声速梯度
    delta_zi =delta_z(size_start:size_end-1);% 各层层厚

    delta_zi(1) = boat_position(i,3) - deepth_a(size_start);% 第一层层厚

    soundspeed_targ_ana = soundspeed(size_start:size_end);% 各层声速
    soundspeed_targ_ana(1) = soundspeed_targ_ana(2) - g_i(1)*delta_zi(1);% 第一层起始位置声速
    size_speed = size(soundspeed_targ_ana,2);

    % 最后一层处理
    delta_zi(size_speed-1) = z_target1 - deepth_a(size_end-1);
    soundspeed_targ_ana(size_speed) = soundspeed_targ_ana(size_speed-1) + g_i(size_speed-1)*delta_zi(size_speed-1);

    sum_t = [];

    for con = 1 : size(angal,2)-2
        % 各层传播所用时间
        simt_i = (asin(p(con) .* (soundspeed_targ_ana(1:size_speed - 1)+g_i .* delta_zi)) ...
                -  asin(p(con) .* soundspeed_targ_ana(1:size_speed - 1)))...
                .*log(1 + g_i .* delta_zi ./ soundspeed_targ_ana(1:size_speed-1))...
                ./(delta_zi .* g_i.^2 * p(con)); 

        % 记录各角度到目标深度传播时间
        sum_t = [sum_t sum(simt_i)]; 
    end
    
    if i == 36
        i = 36;
    end
    for v = 1:5
        if isempty(angal_i)
            angal_c = angal_c;
        else
            angal_c = angal_i;
        end
        
        time_num = find(sum_t > time_delay(i),1);
        
        if abs(sum_t(time_num) - time_delay(i))<delta_tmin
            p_exm = [p_exm sin(angal_c(time_num))./soundspeed_targ(1)];
            angal_i = angal;
            break;
        end

        sum_t = [];
        angal_i = angal_c(time_num-1):(angal_c(time_num)-angal_c(time_num-1))/10:angal_c(time_num);
        
        if isempty(angal_i)
            angal_i = angal_c;
        end

        p_c = sin(angal_i)./soundspeed_targ(1);

        for d = 1:size(angal_i,2)
            
            simt_i1 = (asin(p_c(d) .* (soundspeed_targ_ana(1:size_speed - 1)+g_i .* delta_zi)) ...
                    -  asin(p_c(d) .* soundspeed_targ_ana(1:size_speed - 1)))...
                    .*log(1 + g_i .* delta_zi ./ soundspeed_targ_ana(1:size_speed-1))...
                    ./(delta_zi .* g_i.^2 * p_c(d)); 
            sum_t = [sum_t sum(simt_i1)];
        end
        if v == 5
            time_num = find(sum_t > time_delay(i),1);
            p_exm = [p_exm sin(angal_c(time_num))./soundspeed_targ(1)];
            angal_i = angal;
        end
    end

    % 各层传播水平距离
    simy_i = (sqrt(1 - (p_exm(i) * soundspeed_targ_ana(1:size_speed-1)).^2) ...
            - sqrt(1 - (p_exm(i) * (soundspeed_targ_ana(1:size_speed-1)+g_i.*delta_zi)).^2)) ... 
            ./ (p_exm(i) * g_i + eps);
    %传播水平距离
    sum_y = [sum_y sum(simy_i)];

    distance_meas = sqrt(sum_y.^2 + z_target1^2);
end

error_dis = mean(abs((real_distance - distance_meas'))./real_distance);

distance = distance_meas'.^2; % 平方距离
%% 3. 鲁棒卡尔曼滤波建模
% 状态初始化
noise_level = 10; % 噪声标准差
X = [0; 0; z_target2];     % 初始估计位置(中心区域)
Q = diag([0.01, 0.01, 0.01]);         % 过程噪声(目标静止假设)
R_base = noise_level^2;       % 基础观测噪声
sigma_b_sq = 0;             % 船舶定位误差方差
P = diag([1e6, 1e6, 1e6]); % 初始协方差

% 鲁棒参数配置
k0 = 1.6;                    % Huber阈值系数
max_iter = 1000;               % 最大迭代次数(鲁棒优化)

% 存储变量
estimated_pos = zeros(pos_num, 3);
position_errors = zeros(pos_num, 1);

%% 4. 鲁棒迭代卡尔曼滤波
for i = 1:pos_num
    % 获取当前船舶观测位置
    boat_pos = boat_position(i,:);
    
    % 迭代鲁棒优化
    X_iter = X;
    for iter = 1:max_iter
        % 预测步骤
        P_hat = P + Q;
        
        % 计算雅可比矩阵
        dx = X_iter(1) - boat_pos(1);
        dy = X_iter(2) - boat_pos(2);
        dz = X_iter(3) - boat_pos(3);
        H = [2*dx, 2*dy, 2*dz];
        
        % 计算自适应观测噪声
        dist_sq = dx^2 + dy^2 + dz^2;
        R = R_base + 4*sigma_b_sq*dist_sq;
        
        % 计算残差
        h = dist_sq;
        y = distance(i) - h;
        
        % Huber鲁棒化处理
        S = H*P_hat*H' + R;
        sigma_y = sqrt(S);
        gamma = abs(y)/sigma_y;
        
        if gamma <= k0
            w = 1;
        else
            w = k0/gamma;
        end
        
        % 卡尔曼增益计算
        K = (P_hat*H')/(S);
        
        % 状态更新
        X_new = X_iter + K*(w*y);
        X_new(3) = z_target2;
        % 迭代收敛判断
        if norm(X_new - X_iter) < 1e-3
            break;
        end
        X_iter = X_new;
    end
    
    % 更新协方差
    P = (eye(3) - K*H)*P_hat;
    X = X_iter;
    
    % 记录结果
    estimated_pos(i,:) = X';
    position_errors(i) = norm(X - data_position');
end

%% 高斯牛顿迭代法计算坐标
distance = distance.^0.5;

% 最高迭代次数
n_max = 500;

% 设置门限
delta_max = 10e-9;

% 雅可比行列式J
J = zeros(pos_num-1,3);

% 设定初始值
% position = [0 0 0]';% 参考点坐标轴原点
estimated_pos_real = estimated_pos(pos_num,:)';
position = estimated_pos_real;% 参考点卡尔曼滤波估计点
position_reme = position';
num = 0;
error_num = 0;

for d = 1:100
    distance_pred = zeros(pos_num,1);
    for i = 1:pos_num
        xi = (position(1,1) - boat_position(i,1))^2;
        yi = (position(2,1) - boat_position(i,2))^2;
        zi = (position(3,1) - boat_position(i,3))^2;
        distance_pred(i,1) = (xi + yi + zi)^0.5; % 距离计算
    end
    
    % 计算残差向量r
    r = distance(1:pos_num,1) - distance_pred;
    
    % 计算雅可比行列式
    for i = 1:pos_num
        J(i,1) = position(1,1) - boat_position(i,1);
        J(i,2) = position(2,1) - boat_position(i,2);
        J(i,3) = position(3,1) - boat_position(i,3);
    end
    J = -J./distance_pred;
    
    % 计算误差
    delta = -((J'*J)^-1)*(J'*r); 
    absdelta = norm(delta);

    % 比较门限
    if absdelta < delta_max
        break;
    else
        position = position + delta;
        delta(3) = 0;
        position(3) = z_target2;
        position_reme = [position_reme;position'];
    end

    if d == 100
        error_num = 1;
        break;
    end
end
toc
fprintf('真实坐标: (%.2f, %.2f, %.2f)\n', data_position);
fprintf('初始估计: (%.2f, %.2f, %.2f) 误差: %.2f米\n', ...
    estimated_pos_real, norm(estimated_pos_real - data_position'));
fprintf('最终估计: (%.2f, %.2f, %.2f) 误差: %.2f米\n', ...
    position, norm(position - data_position'));
fprintf('测距误差: %.2f %% ', error_dis*100);
%% 4. 结果可视化
figure;
% 绘制船的位置
boat_pos_disp = boat_position;

% boat_posreal_disp = [boat_position_real;boat_position_real(1,:)];
% scatter3(boat_pos_disp(:,1), boat_pos_disp(:,2), boat_pos_disp(:,3),'.','MarkerFaceColor','b');
plot3(boat_pos_disp(:,1), boat_pos_disp(:,2), boat_pos_disp(:,3),'.-','MarkerFaceColor','b');
hold on;
% plot3(position_reme(:,1), position_reme(:,2), position_reme(:,3),'.--','MarkerFaceColor','y');
% legend('船舶测量位置','高斯牛顿迭代过程');
% scatter3(boat_posreal_disp(:,1), boat_posreal_disp(:,2), boat_posreal_disp(:,3), ...
%          100, 'r', 'x', 'LineWidth',1, 'DisplayName', '船舶真实位置');
% 绘制真实目标位置
scatter3(data_position(1), data_position(2), data_position(3), ...
         100, 'g', 'pentagram', 'filled', 'DisplayName', '真实位置');
% 绘制估计位置
scatter3(position(1), position(2), position(3), ...
         100, 'r', 'x', 'LineWidth',2, 'DisplayName', '估计位置');
title("等梯度声线跟踪");
grid on % 显示主网格线

%% 4. 结果可视化
figure;
% 绘制船的位置
boat_pos_disp = boat_position;
% boat_posreal_disp = [boat_position_real;boat_position_real(1,:)];
% scatter3(boat_pos_disp(:,1), boat_pos_disp(:,2), boat_pos_disp(:,3),'.','MarkerFaceColor','b');
plot(boat_pos_disp(:,1), boat_pos_disp(:,2), '.-','MarkerFaceColor','b');
hold on;
% plot3(position_reme(:,1), position_reme(:,2), position_reme(:,3),'.--','MarkerFaceColor','y');
% legend('船舶测量位置','高斯牛顿迭代过程');
% scatter3(boat_posreal_disp(:,1), boat_posreal_disp(:,2), boat_posreal_disp(:,3), ...
%          100, 'r', 'x', 'LineWidth',1, 'DisplayName', '船舶真实位置');
% 绘制真实目标位置
scatter(data_position(1), data_position(2), ...
         100, 'r', '^', 'filled', 'DisplayName', '真实位置');
% 绘制估计位置
% scatter3(position(1), position(2), position(3), ...
%          100, 'r', 'x', 'LineWidth',2, 'DisplayName', '估计位置');
title("等梯度声线跟踪");
grid on % 显示主网格线