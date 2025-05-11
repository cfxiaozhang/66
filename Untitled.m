%% 改进蚁群算法求解带时间窗的车辆路径问题 (VRPTW)
clear;
clc;
close all;

%% 参数设置
n = 50;                 % 客户数量
m = 20;                 % 车辆数量
capacity = 100;          % 车辆容量
depot = [0, 0];         % 仓库坐标
F = 300;                 % 车辆固定使用成本
C = 2;                  % 单位距离运输成本
mu = 0.001;               % 单位时间损耗系数
P = 10000;                 % 单位货物价值
G = 0.5;                  % 冷链设备单位时间成本
late_penalty = 2;       % 迟到惩罚系数
speed = 1;              % 车辆速度

%% 随机生成客户数据
rng(0);
customers.x = rand(n, 1) * 100;
customers.y = rand(n, 1) * 100;
customers.demand = randi([1, 10], n, 1);
customers.ready_time = randi([0, 50], n, 1);
customers.due_time = customers.ready_time + randi([20, 50], n, 1);
customers.service_time = randi([5, 10], n, 1);

%% 蚁群算法参数
num_ants = 100;
max_iterations = 300;
alpha0 = 1;     % 初始信息素重要程度
beta0 = 5;      % 初始启发因子重要程度
rho0 = 0.8;     % 初始信息素挥发系数
Q = 100;
elitist_weight = 10;  % 精英蚂蚁加权因子

tau_max_init = 1;     % 初始最大信息素
tau_min_init = 1e-4;  % 初始最小信息素

alpha = alpha0;
beta = beta0;
rho = rho0;
tau_max = tau_max_init;
tau_min = tau_min_init;

%% 距离矩阵
distance = zeros(n+1, n+1);
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
pheromone = tau_max * ones(n+1, n+1);
for i = 1:n+1
    pheromone(i, i) = 0;
end

%% 记录变量
best_cost = inf;
best_routes = [];
best_iteration = 0;
iter_best_cost = zeros(max_iterations, 1);

%% 主循环
for iter = 1:max_iterations
    % *** 自适应参数调整 ***
    alpha = alpha0 + (3 - alpha0) * iter / max_iterations;  % alpha渐增
    beta  = beta0  - (beta0 - 1) * iter / max_iterations;   % beta渐减
    rho   = rho0   * (1 - iter / (2*max_iterations));        % 逐步减小rho
    
    % *** 最大最小信息素限制 ***
    if iter == 1
        tau_max = tau_max_init;
        tau_min = tau_min_init;
    else
        tau_max = 1 / (rho * best_cost);
        tau_min = tau_max / (2 * n);
    end

    all_routes = cell(num_ants, 1);
    all_costs = zeros(num_ants, 1);
    
    % 蚂蚁逐个搜索
    for ant = 1:num_ants
        unvisited = 2:n+1;
        routes = cell(m, 1);
        vehicle_load = zeros(m, 1);
        vehicle_time = zeros(m, 1);
        route_index = ones(m, 1);
        for k = 1:m
            current_node = 1;
            routes{k}(route_index(k)) = current_node;
            while ~isempty(unvisited) && vehicle_load(k) < capacity
                probabilities = zeros(size(unvisited));
                for i = 1:length(unvisited)
                    j = unvisited(i);
                    if vehicle_load(k) + customers.demand(j-1) <= capacity
                        travel_time = distance(current_node, j) / speed;
                        arrival_time = vehicle_time(k) + travel_time;
                        wait_time = max(0, customers.ready_time(j-1) - arrival_time);
                        departure_time = arrival_time + wait_time + customers.service_time(j-1);
                        late_penalty_cost = late_penalty * max(0, departure_time - customers.due_time(j-1));
                        eta = 1 / (distance(current_node, j) + late_penalty_cost + 0.001);
                        probabilities(i) = pheromone(current_node, j)^alpha * eta^beta;
                    end
                end
                if sum(probabilities) == 0
                    break;
                end
                probabilities = probabilities / sum(probabilities);
                r = rand;
                cumulative = 0;
                next_index = 1;
                for i = 1:length(probabilities)
                    cumulative = cumulative + probabilities(i);
                    if r <= cumulative
                        next_index = i; break;
                    end
                end
                next_node = unvisited(next_index);
                vehicle_load(k) = vehicle_load(k) + customers.demand(next_node-1);
                travel_time = distance(current_node, next_node) / speed;
                arrival_time = vehicle_time(k) + travel_time;
                wait_time = max(0, customers.ready_time(next_node-1) - arrival_time);
                vehicle_time(k) = arrival_time + wait_time + customers.service_time(next_node-1);
                route_index(k) = route_index(k) + 1;
                routes{k}(route_index(k)) = next_node;
                unvisited(next_index) = [];
                current_node = next_node;
            end
            route_index(k) = route_index(k) + 1;
            routes{k}(route_index(k)) = 1;
        end
        
        % *** 两两交换法(2-opt)局部搜索 ***
        for k = 1:m
            route = routes{k};
            if length(route) > 3
                routes{k} = two_opt(route, distance);
            end
        end
        
        % 计算总代价
        total_cost = 0;
        for k = 1:m
            route = routes{k};
            if length(route) > 2
                vehicle_time = 0;
                vehicle_distance = 0;
                damage_cost = 0;
                equipment_cost = 0;
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
                        late_cost = late_penalty * max(0, departure_time - customers.due_time(j));
                        total_cost = total_cost + late_cost;
                        vehicle_time = departure_time;
                    else
                        vehicle_time = vehicle_time + travel_time;
                    end
                end
                total_cost = total_cost + vehicle_distance * C + F;
                total_cost = total_cost + damage_cost + equipment_cost;
            end
        end
        
        all_routes{ant} = routes;
        all_costs(ant) = total_cost;
        if total_cost < best_cost
            best_cost = total_cost;
            best_routes = routes;
            best_iteration = iter;
        end
    end
    
    % 记录本轮最优
    [min_cost, idx] = min(all_costs);
    iter_best_cost(iter) = min_cost;

    % 信息素挥发
    pheromone = (1 - rho) * pheromone;
    
    % 蚂蚁路径更新
    for ant = 1:num_ants
        routes = all_routes{ant};
        for k = 1:m
            route = routes{k};
            for i = 1:length(route)-1
                pheromone(route(i), route(i+1)) = pheromone(route(i), route(i+1)) + Q / all_costs(ant);
            end
        end
    end
    % *** 精英蚂蚁系统：对历史最优路径额外强化信息素 ***
    for k = 1:m
        route = best_routes{k};
        for i = 1:length(route)-1
            pheromone(route(i), route(i+1)) = pheromone(route(i), route(i+1)) + elitist_weight * Q / best_cost;
        end
    end
    % *** 最大-最小信息素限制（MMAS） ***
    pheromone(pheromone > tau_max) = tau_max;
    pheromone(pheromone < tau_min) = tau_min;

    % 显示进度
    if mod(iter, 10) == 0
        fprintf('迭代次数: %d, 当前最优总代价: %.2f\n', iter, best_cost);
    end
end

%% 输出结果
fprintf('最优解出现在第 %d 代。\n', best_iteration);
fprintf('总代价: %.2f\n', best_cost);

fprintf('路径方案:\n');
for k = 1:m
    route = best_routes{k};
    if length(route) > 2
        fprintf('车辆 %d: ', k);
        for i = 1:length(route)
            if route(i) == 1
                fprintf('仓库');
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

% 各项成本统计
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
fprintf('运输成本(距离*C): %.2f\n', total_distance * C);
fprintf('车辆固定成本: %.2f\n', total_fixed_cost);
fprintf('迟到惩罚: %.2f\n', total_late_penalty);
fprintf('损耗成本: %.2f\n', total_damage_cost);
fprintf('冷链设备成本: %.2f\n', total_equipment_cost);
fprintf('总代价: %.2f\n', total_distance * C + total_fixed_cost + total_late_penalty + total_damage_cost + total_equipment_cost);

%% 收敛曲线
figure;
plot(1:max_iterations, iter_best_cost, 'LineWidth', 2);
title('蚁群算法收敛曲线');
xlabel('迭代次数');
ylabel('总代价');
grid on;

%% 路径可视化
figure; hold on;
plot(depot(1), depot(2), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
text(depot(1)+1, depot(2)+1, '仓库');
for i = 1:n
    plot(customers.x(i), customers.y(i), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    text(customers.x(i)+1, customers.y(i)+1, num2str(i));
end
colors = lines(m);
for k = 1:m
    route = best_routes{k};
    if length(route) > 2
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
title('最优车辆路径');
xlabel('X坐标');
ylabel('Y坐标');
grid on;
axis equal;
legend('仓库', '客户', '车辆路径');
hold off;