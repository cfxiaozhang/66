%% 蚁群算法求解带软时间窗的车辆路径问题(VRPTW)
clear;
clc;
close all;

%% 问题参数设置
n = 50;                 % 客户数量
m = 20;                  % 车辆数量
capacity = 50;          % 车辆容量
depot = [0, 0];         % 配送中心坐标
F = 50;                 % 车辆固定使用成本（文档中的F）
C = 2;                  % 单位距离行驶成本（文档中的C）
mu = 0.1;               % 单位时间货损系数（文档中的μ）
P = 10;                 % 单位卷烟价值（文档中的P）
G = 5;                  % 温湿度设备单位时间成本（文档中的G）
late_penalty = 2;       % 晚到惩罚系数（文档中仅晚到有惩罚）
speed = 1;              % 车辆行驶速度

%% 生成随机客户数据
rng(0); % 设置随机种子，保证结果可重现
customers.x = rand(n, 1) * 100;      % 客户x坐标
customers.y = rand(n, 1) * 100;      % 客户y坐标
customers.demand = randi([1, 10], n, 1); % 客户需求
customers.ready_time = randi([0, 50], n, 1); % 客户最早服务时间
customers.due_time = customers.ready_time + randi([20, 50], n, 1); % 客户最晚服务时间
customers.service_time = randi([5, 10], n, 1); % 客户服务时间

%% 蚁群算法参数设置
num_ants = 100;          % 蚂蚁数量
max_iterations = 300;   % 最大迭代次数
alpha = 1;              % 信息素重要程度因子
beta = 2;               % 启发式因子
rho = 0.5;              % 信息素挥发系数
Q = 100;                % 信息素更新强度

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

%% 初始化信息素矩阵
pheromone = ones(n+1, n+1); % 初始信息素
for i = 1:n+1
    pheromone(i, i) = 0; % 自身到自身的信息素为0
end

%% 初始化最佳解
best_cost = inf;
best_routes = [];
best_iteration = 0;
iter_best_cost = zeros(max_iterations, 1); % 记录每次迭代的最优解

%% 蚁群算法主循环
for iter = 1:max_iterations
    % 初始化所有蚂蚁的路径和成本
    all_routes = cell(num_ants, 1);
    all_costs = zeros(num_ants, 1);
    
    % 每只蚂蚁构建路径
    for ant = 1:num_ants
        % 初始化
        unvisited = 2:n+1; % 未访问的客户（2到n+1，1是配送中心）
        routes = cell(m, 1); % m辆车的路径
        vehicle_load = zeros(m, 1); % 每辆车的当前载重
        vehicle_time = zeros(m, 1); % 每辆车的当前时间
        route_index = ones(m, 1); % 每辆车当前路径的索引
        
        % 为每辆车分配路径
        for k = 1:m
            current_node = 1; % 从配送中心开始
            routes{k}(route_index(k)) = current_node;
            
            % 当还有未访问的客户且车辆未满载时继续
            while ~isempty(unvisited) && vehicle_load(k) < capacity
                % 计算转移概率
                probabilities = zeros(size(unvisited));
                for i = 1:length(unvisited)
                    j = unvisited(i);
                    if vehicle_load(k) + customers.demand(j-1) <= capacity
                        % 计算时间相关参数
                        travel_time = distance(current_node, j) / speed;
                        arrival_time = vehicle_time(k) + travel_time;
                        wait_time = max(0, customers.ready_time(j-1) - arrival_time);
                        departure_time = arrival_time + wait_time + customers.service_time(j-1);
                        late_penalty_cost = late_penalty * max(0, departure_time - customers.due_time(j-1));
                        
                        % 计算启发式信息
                        eta = 1 / (distance(current_node, j) + late_penalty_cost + 0.001);
                        
                        % 计算转移概率
                        probabilities(i) = pheromone(current_node, j)^alpha * eta^beta;
                    end
                end
                
                % 如果没有可行的客户，结束当前车辆的路径
                if sum(probabilities) == 0
                    break;
                end
                
                % 归一化概率
                probabilities = probabilities / sum(probabilities);
                
                % 轮盘赌选择下一个客户
                r = rand;
                cumulative = 0;
                next_index = 1;
                for i = 1:length(probabilities)
                    cumulative = cumulative + probabilities(i);
                    if r <= cumulative
                        next_index = i;
                        break;
                    end
                end
                next_node = unvisited(next_index);
                
                % 更新车辆状态
                vehicle_load(k) = vehicle_load(k) + customers.demand(next_node-1);
                travel_time = distance(current_node, next_node) / speed;
                arrival_time = vehicle_time(k) + travel_time;
                wait_time = max(0, customers.ready_time(next_node-1) - arrival_time);
                vehicle_time(k) = arrival_time + wait_time + customers.service_time(next_node-1);
                
                % 更新路径
                route_index(k) = route_index(k) + 1;
                routes{k}(route_index(k)) = next_node;
                
                % 从未访问客户中移除
                unvisited(next_index) = [];
                
                % 更新当前节点
                current_node = next_node;
            end
            
            % 车辆返回配送中心
            route_index(k) = route_index(k) + 1;
            routes{k}(route_index(k)) = 1;
        end
        
        % 计算总成本（使用修改后的目标函数）
        total_cost = 0;
        vehicle_used = 0;       % 记录使用的车辆数
        for k = 1:m
            route = routes{k};
            if length(route) > 2 % 至少有配送中心和一个客户
                vehicle_used = vehicle_used + 1;
                vehicle_time = 0;
                vehicle_distance = 0;
                damage_cost = 0;
                equipment_cost = 0;
                
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
                        wait_time = max(0, customers.ready_time(j) - arrival_time);
                        departure_time = arrival_time + wait_time + customers.service_time(j);
                        late_cost = late_penalty * max(0, departure_time - customers.due_time(j));
                        total_cost = total_cost + late_cost;
                        
                        vehicle_time = departure_time;
                    else
                        vehicle_time = vehicle_time + travel_time; % 返回配送中心的时间
                    end
                end
                
                % 累加固定成本和行驶成本
                total_cost = total_cost + vehicle_distance * C + F;
                total_cost = total_cost + damage_cost + equipment_cost;
            end
        end
        
        % 保存当前蚂蚁的路径和成本
        all_routes{ant} = routes;
        all_costs(ant) = total_cost;
        
        % 更新全局最优解
        if total_cost < best_cost
            best_cost = total_cost;
            best_routes = routes;
            best_iteration = iter;
        end
    end
    
    % 记录当前迭代的最优解
    [min_cost, ~] = min(all_costs);
    iter_best_cost(iter) = min_cost;
    
    % 信息素更新
    pheromone = (1 - rho) * pheromone; % 信息素挥发
    
    % 每只蚂蚁更新信息素
    for ant = 1:num_ants
        routes = all_routes{ant};
        for k = 1:m
            route = routes{k};
            for i = 1:length(route)-1
                pheromone(route(i), route(i+1)) = pheromone(route(i), route(i+1)) + Q / all_costs(ant);
            end
        end
    end
    
    % 显示进度
    if mod(iter, 10) == 0
        fprintf('迭代次数: %d, 最优成本: %.2f\n', iter, best_cost);
    end
end

%% 输出结果
fprintf('最优解在第 %d 次迭代找到\n', best_iteration);
fprintf('总成本: %.2f\n', best_cost);

% 显示最优路径
fprintf('最优路径:\n');
for k = 1:m
    route = best_routes{k};
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
total_distance = 0;
total_fixed_cost = 0;
total_late_penalty = 0;
total_damage_cost = 0;
total_equipment_cost = 0;

for k = 1:m
    route = best_routes{k};
    if length(route) > 2
        vehicle_time = 0;
        vehicle_distance = 0;
        damage_cost = 0;
        equipment_cost = 0;
        late_cost = 0;
        
        for i = 1:length(route)-1
            current_node = route(i);
            next_node = route(i+1);
            
            vehicle_distance = vehicle_distance + distance(current_node, next_node);
            travel_time = distance(current_node, next_node) / speed;
            equipment_cost = equipment_cost + travel_time * G;
            
            if next_node ~= 1
                j = next_node - 1;
                damage_cost = damage_cost + travel_time * mu * P * customers.demand(j);
                
                vehicle_time = vehicle_time + travel_time;
                arrival_time = vehicle_time;
                wait_time = max(0, customers.ready_time(j) - arrival_time);
                departure_time = arrival_time + wait_time + customers.service_time(j);
                late_cost = late_cost + late_penalty * max(0, departure_time - customers.due_time(j));
                vehicle_time = departure_time;
            else
                vehicle_time = vehicle_time + travel_time;
            end
        end
        
        total_distance = total_distance + vehicle_distance;
        total_fixed_cost = total_fixed_cost + F;
        total_late_penalty = total_late_penalty + late_cost;
        total_damage_cost = total_damage_cost + damage_cost;
        total_equipment_cost = total_equipment_cost + equipment_cost;
    end
end

fprintf('行驶成本 (距离×C): %.2f\n', total_distance * C);
fprintf('车辆固定成本: %.2f\n', total_fixed_cost);
fprintf('晚到惩罚成本: %.2f\n', total_late_penalty);
fprintf('货损成本: %.2f\n', total_damage_cost);
fprintf('温湿度设备成本: %.2f\n', total_equipment_cost);
fprintf('总成本: %.2f\n', total_distance * C + total_fixed_cost + total_late_penalty + total_damage_cost + total_equipment_cost);

%% 可视化
% 迭代过程图
figure;
plot(1:max_iterations, iter_best_cost, 'LineWidth', 2);
title('蚁群算法迭代过程');
xlabel('迭代次数');
ylabel('总成本');
grid on;

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
    route = best_routes{k};
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

title('最优车辆路径方案');
xlabel('X坐标');
ylabel('Y坐标');
grid on;
axis equal;
legend('配送中心', '客户', '车辆路径');
hold off;    