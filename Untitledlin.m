%% 最近邻算法求解带软时间窗的车辆路径问题(VRPTW)
clear;
clc;
close all;

%% 问题参数设置
n = 50;                 % 客户数量
m = 25;                  % 车辆数量
capacity = 200;          % 车辆容量
depot = [35, 35];         % 配送中心坐标
F = 300;                 % 车辆固定使用成本
C = 2;                  % 单位距离行驶成本
mu = 0.001;               % 单位时间货损系数
P = 10000;                 % 单位卷烟价值
G = 0.5;                  % 温湿度设备单位时间成本
late_penalty = 2;       % 晚到惩罚系数
speed = 60;              % 车辆行驶速度

%% 生成随机客户数据
rng(0); % 设置随机种子，保证结果可重现
customers.x = xlsread('客户点.xlsx','B11:B60');      % 客户x坐标
customers.y = xlsread('客户点.xlsx','C11:C60');      % 客户y坐标
customers.demand = xlsread('客户点.xlsx','D11:D60'); % 客户需求
customers.ready_time = xlsread('客户点.xlsx','E11:E60'); % 客户最早服务时间
customers.due_time = xlsread('客户点.xlsx','F11:F60'); % 客户最晚服务时间
customers.service_time = xlsread('客户点.xlsx','G11:G60'); % 客户服务时间

%% 计算距离矩阵
distance = zeros(n+1, n+1); % 包括配送中心
for i = 1:n+1
    for j = 1:n+1
        if i ~= j
            if i == 1
                x1 = depot(1);
                y1 = depot(2);
            else
                x1 = customers.x(i-1);
                y1 = customers.y(i-1);
            end
            
            if j == 1
                x2 = depot(1);
                y2 = depot(2);
            else
                x2 = customers.x(j-1);
                y2 = customers.y(j-1);
            end
            
            distance(i, j) = sqrt((x1 - x2)^2 + (y1 - y2)^2);
        end
    end
end

%% 最近邻算法求解
unvisited = 2:n+1; % 未访问的客户（2到n+1，1是配送中心）
routes = cell(m, 1); % m辆车的路径
vehicle_load = zeros(m, 1); % 每辆车的当前载重
vehicle_time = zeros(m, 1); % 每辆车的当前时间
total_cost = 0; % 总成本

for k = 1:m
    current_node = 1; % 从配送中心开始
    routes{k} = [current_node];
    
    while ~isempty(unvisited)
        nearest_node = 0;
        min_distance = inf;
        
        % 寻找最近的可行客户
        for i = 1:length(unvisited)
            j = unvisited(i);
            
            % 检查容量约束
            if vehicle_load(k) + customers.demand(j-1) > capacity
                continue;
            end
            
            % 计算到达时间和时间窗约束
            travel_time = distance(current_node, j) / speed;
            arrival_time = vehicle_time(k) + travel_time;
            
            % 计算晚到惩罚（软时间窗）
            late_penalty_cost = late_penalty * max(0, arrival_time - customers.due_time(j-1));
            
            % 计算总成本（距离+晚到惩罚）
            total_cost_candidate = distance(current_node, j) + late_penalty_cost;
            
            % 如果更近且满足约束，更新最近节点
            if total_cost_candidate < min_distance
                min_distance = total_cost_candidate;
                nearest_node = j;
            end
        end
        
        % 如果没有找到可行客户，结束当前车辆的路径
        if nearest_node == 0
            break;
        end
        
        % 更新车辆状态
        vehicle_load(k) = vehicle_load(k) + customers.demand(nearest_node-1);
        travel_time = distance(current_node, nearest_node) / speed;
        vehicle_time(k) = vehicle_time(k) + travel_time;
        wait_time = max(0, customers.ready_time(nearest_node-1) - vehicle_time(k));
        vehicle_time(k) = vehicle_time(k) + wait_time + customers.service_time(nearest_node-1);
        
        % 更新路径
        routes{k} = [routes{k}, nearest_node];
        
        % 从未访问客户中移除
        unvisited(find(unvisited == nearest_node)) = [];
        
        % 更新当前节点
        current_node = nearest_node;
    end
    
    % 车辆返回配送中心
    if length(routes{k}) > 1 % 如果车辆有服务客户
        routes{k} = [routes{k}, 1];
        
        % 计算车辆返回配送中心的时间
        travel_time = distance(current_node, 1) / speed;
        vehicle_time(k) = vehicle_time(k) + travel_time;
    end
end

%% 计算总成本
total_cost = 0;
total_distance = 0;
total_fixed_cost = 0;
total_late_penalty = 0;
total_damage_cost = 0;
total_equipment_cost = 0;

for k = 1:m
    route = routes{k};
    if length(route) > 2 % 至少有配送中心和一个客户
        vehicle_time = 0;
        vehicle_distance = 0;
        damage_cost = 0;
        equipment_cost = 0;
        late_cost = 0;
        
        total_fixed_cost = total_fixed_cost + F;
        
        for i = 1:length(route)-1
            current_node = route(i);
            next_node = route(i+1);
            
            % 行驶成本
            vehicle_distance = vehicle_distance + distance(current_node, next_node);
            
            % 运输时间（用于货损和设备成本）
            travel_time = distance(current_node, next_node) / speed;
            equipment_cost = equipment_cost + travel_time * G;
            
            if next_node ~= 1 % 非配送中心节点（客户节点）
                j = next_node - 1; % 客户索引
                damage_cost = damage_cost + travel_time * mu * P * customers.demand(j);
                
                % 时间窗惩罚（仅晚到惩罚）
                vehicle_time = vehicle_time + travel_time;
                arrival_time = vehicle_time;
                late_cost = late_cost + late_penalty * max(0, arrival_time - customers.due_time(j));
                
                % 更新车辆时间
                wait_time = max(0, customers.ready_time(j) - arrival_time);
                vehicle_time = arrival_time + wait_time + customers.service_time(j);
            else
                vehicle_time = vehicle_time + travel_time; % 返回配送中心的时间
            end
        end
        
        total_distance = total_distance + vehicle_distance;
        total_late_penalty = total_late_penalty + late_cost;
        total_damage_cost = total_damage_cost + damage_cost;
        total_equipment_cost = total_equipment_cost + equipment_cost;
    end
end

% 计算总成本
total_cost = total_distance * C + total_fixed_cost + total_late_penalty + total_damage_cost + total_equipment_cost;

%% 输出结果
fprintf('最近邻算法求解结果:\n');
fprintf('总成本: %.2f\n', total_cost);

% 显示最优路径
fprintf('最优路径:\n');
for k = 1:m
    route = routes{k};
    if length(route) > 2 % 至少有配送中心和一个客户
        fprintf('车辆 %d: ', k);
        for i = 1:length(route)
            if route(i) == 1
                fprintf('配送中心');
            else
                fprintf('客户 %d', route(i)-1);
            end
            if i < length(route)
                fprintf(' -> ');
            end
        end
        fprintf('\n');
    end
end

% 分解显示各项成本
fprintf('\n成本分解:\n');
fprintf('行驶成本 (距离×C): %.2f\n', total_distance * C);
fprintf('车辆固定成本: %.2f\n', total_fixed_cost);
fprintf('晚到惩罚成本: %.2f\n', total_late_penalty);
fprintf('货损成本: %.2f\n', total_damage_cost);
fprintf('温湿度设备成本: %.2f\n', total_equipment_cost);
fprintf('总成本: %.2f\n', total_cost);

%% 可视化
% 最优路径图
figure;
hold on;
% 绘制配送中心
plot(depot(1), depot(2), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
text(depot(1)+1, depot(2)+1, '配送中心');

% 绘制客户点
for i = 1:n
    plot(customers.x(i), customers.y(i), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    text(customers.x(i)+1, customers.y(i)+1, num2str(i));
end

% 绘制路径
colors = lines(m);
for k = 1:m
    route = routes{k};
    if length(route) > 2 % 至少有配送中心和一个客户
        x = zeros(length(route), 1);
        y = zeros(length(route), 1);
        for i = 1:length(route)
            if route(i) == 1
                x(i) = depot(1);
                y(i) = depot(2);
            else
                x(i) = customers.x(route(i)-1);
                y(i) = customers.y(route(i)-1);
            end
        end
        plot(x, y, 'Color', colors(k,:), 'LineWidth', 1.5);
        text(mean(x), mean(y), ['车辆 ', num2str(k)], 'Color', colors(k,:));
    end
end

title('最近邻算法最优车辆路径方案');
xlabel('X坐标');
ylabel('Y坐标');
grid on;
axis equal;
legend('配送中心', '客户', '车辆路径');
hold off;    