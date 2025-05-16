clear; clc; close all; % 清理工作区，命令窗口，关闭所有图形

%% 1. 问题定义
% 客户数据: x坐标, y坐标, 需求量, 最早到达时间, 最晚到达时间, 服务时长, 货物单位价值, 货损率
s = rng(1); % 设置随机数种子以便结果可复现

num_customers = 50; % 客户数量
customer_data = struct('x', num2cell(xlsread('客户点.xlsx','A2:AX2')), ... % x坐标
                       'y', num2cell(xlsread('客户点.xlsx','A3:AX3')), ... % y坐标
                       'demand', num2cell(xlsread('客户点.xlsx','A4:AX4')), ... % 需求量
                       'early_tw', num2cell(xlsread('客户点.xlsx','A5:AX5')), ... % 最早到达时间
                       'late_tw', num2cell(xlsread('客户点.xlsx','A6:AX6')), ...% 最晚到达时间
                       'service_time', num2cell(xlsread('客户点.xlsx','A7:AX7')), ... % 服务时长
                       'cargo_value', num2cell(1000 * ones(1, num_customers)), ...    % 货物单位价值
                       'damage_rate', num2cell(0.001 * ones(1, num_customers)));   % 货损率 (例如1%)

% 确保最晚时间窗晚于 (最早时间窗 + 服务时长)，尽管软时间窗会处理，但这样设置更合理
for i=1:num_customers
    customer_data(i).late_tw = max(customer_data(i).late_tw, customer_data(i).early_tw + customer_data(i).service_time + 10); % 加10是为了留些余地
end

% 配送中心 (Depot) 数据
depot_data.x = 35;         % Depot x坐标
depot_data.y = 35;         % Depot y坐标
depot_data.early_tw = 0;   % Depot 最早可出发时间
depot_data.late_tw = 230;  % Depot 最晚运营结束时间 (车辆必须在此之前返回)

% 参数设置
params.vehicle_capacity = 100;         % 车辆容量
params.fixed_vehicle_cost = 300;      % 每辆车的固定成本
params.cost_per_distance = 3;         % 单位距离的行驶成本
params.penalty_per_late_time = 3;     % 单位晚到时间的惩罚成本
params.unit_refrigeration_cost = 1; % 单位制冷时间的成本

% 整合Depot和客户数据，方便计算距离矩阵和索引
% 节点1: Depot, 节点 2 到 N+1: 客户
all_nodes_coords = [depot_data.x, depot_data.y; [customer_data.x]', [customer_data.y]']; % 所有节点坐标
all_nodes_demand = [0; [customer_data.demand]']; % 所有节点需求 (Depot为0)
all_nodes_early_tw = [depot_data.early_tw; [customer_data.early_tw]']; % 所有节点最早时间窗
all_nodes_late_tw = [depot_data.late_tw; [customer_data.late_tw]'];   % 所有节点最晚时间窗
all_nodes_service_time = [0; [customer_data.service_time]']; % 所有节点服务时长 (Depot为0)

num_total_nodes = size(all_nodes_coords, 1); % 总节点数

% 计算距离矩阵
dist_matrix = zeros(num_total_nodes, num_total_nodes);
for i = 1:num_total_nodes
    for j = i+1:num_total_nodes
        dist_matrix(i,j) = norm(all_nodes_coords(i,:) - all_nodes_coords(j,:));
        dist_matrix(j,i) = dist_matrix(i,j);
    end
end

%% 2. 蚁群算法 (ACO) 参数
num_ants = 100;        % 蚂蚁数量
num_iterations = 300; % 迭代次数
alpha = 1;            % 信息素重要程度因子
beta = 3;             % 启发式信息 (距离) 重要程度因子
rho = 0.5;            % 信息素挥发率
Q = 100;              % 信息素增加强度系数

% 启发式信息 (eta = 1/距离)
eta = 1 ./ dist_matrix;
eta(isinf(eta) | isnan(eta)) = 0; % 处理除以零的情况 (同一个节点距离为0)

% 初始化信息素矩阵 (tau)
tau = ones(num_total_nodes, num_total_nodes) * 1e-3; % 初始信息素设为一个较小的值

best_solution_routes = {};         % 存储全局最优解的路径
best_solution_cost = inf;          % 存储全局最优解的成本
best_costs_history = zeros(1, num_iterations); % 记录每次迭代的最优成本

%% 3. ACO 主循环
fprintf('CVRPTW的蚁群算法开始...\n');
for iter = 1:num_iterations
    fprintf('迭代 %d/%d\n', iter, num_iterations);

    ant_solutions = cell(num_ants, 1); % 存储当前迭代中每只蚂蚁的解 (路径集合)
    ant_costs = zeros(num_ants, 1);    % 存储当前迭代中每只蚂蚁解的成本

    for ant = 1:num_ants % 对于每只蚂蚁
        ant_routes = {}; % 当前蚂蚁构建的路径集合
        visited_customers = false(1, num_customers); % 标记客户是否已被访问 (客户索引从1到num_customers, 对应节点索引2到N+1)
        num_vehicles_ant = 0; % 当前蚂蚁使用的车辆数

        while ~all(visited_customers) % 当仍有客户未被服务时
            num_vehicles_ant = num_vehicles_ant + 1;
            current_route = [1]; % 新路径从Depot (节点1) 开始
            current_load = 0;    % 当前车辆载重
            current_time = all_nodes_early_tw(1); % 车辆从Depot出发的时间
            current_node = 1;    % 当前所在节点

            departure_time_from_node = current_time; % 从当前节点出发的时间

            while true % 构建单条车辆路径
                possible_next_nodes = []; % 可能的下一个客户节点列表
                probabilities = [];       % 对应节点的选择概率

                % 寻找下一个可服务的客户节点
                for next_node_idx = 2:num_total_nodes % 遍历所有客户节点 (索引2到N+1)
                    customer_actual_idx = next_node_idx - 1; % 对应 customer_data 中的索引
                    if ~visited_customers(customer_actual_idx) && ... % 如果该客户未被全局访问过
                       (current_load + all_nodes_demand(next_node_idx) <= params.vehicle_capacity) % 且不超过车辆容量

                        % 时间窗可行性初步检查 (更详细的检查在成本计算函数中)
                        travel_time_to_next = dist_matrix(current_node, next_node_idx);
                        arrival_at_next = departure_time_from_node + travel_time_to_next;

                        % 在下一个节点的最早服务开始时间
                        service_start_at_next = max(arrival_at_next, all_nodes_early_tw(next_node_idx));

                        % 检查从下一个节点服务完毕后返回depot是否会超过depot的最晚运营时间
                        departure_from_next = service_start_at_next + all_nodes_service_time(next_node_idx);
                        time_to_return_depot = dist_matrix(next_node_idx, 1); % 从next_node_idx返回depot的时间

                        if departure_from_next + time_to_return_depot <= all_nodes_late_tw(1) % Depot的最晚关闭时间
                             possible_next_nodes = [possible_next_nodes, next_node_idx];
                        end
                    end
                end

                if isempty(possible_next_nodes) % 如果没有可选择的下一个节点
                    break; % 结束当前车辆的路径构建 (可能车辆已满，或剩余客户不满足时间约束)
                end

                % 计算选择下一个节点的概率
                prob_numerator = []; % 存储概率公式的分子部分
                for nn_idx = 1:length(possible_next_nodes)
                    node_cand = possible_next_nodes(nn_idx); % 候选节点
                    prob_numerator(nn_idx) = (tau(current_node, node_cand)^alpha) * ... % 信息素浓度
                                             (eta(current_node, node_cand)^beta);      % 启发式信息
                end

                if sum(prob_numerator) == 0 % 避免除以零 (如果所有分子都为0)
                    % 这种情况可能发生在eta为0 (例如beta很高且距离很大) 或tau极小。
                    % 简单处理：随机选择一个，或中止当前车辆路径。
                     if isempty(possible_next_nodes)
                        break;
                     else % 如果有多个概率为0的选项，随机选一个
                        selected_idx_in_possible = randi(length(possible_next_nodes));
                     end
                else
                    probabilities = prob_numerator / sum(prob_numerator); % 归一化得到概率
                    % 轮盘赌选择下一个节点
                    cumulative_prob = cumsum(probabilities);
                    r_select = rand();
                    selected_idx_in_possible = find(r_select <= cumulative_prob, 1, 'first');
                end

                if isempty(selected_idx_in_possible) % 如果没有选出节点 (理论上不应发生，除非possible_next_nodes为空或sum(prob_num)<=0)
                     break; % 中止当前车辆路径
                end

                next_chosen_node = possible_next_nodes(selected_idx_in_possible); % 选中的下一个节点

                % 更新路径、载重、时间等信息
                current_route = [current_route, next_chosen_node];
                current_load = current_load + all_nodes_demand(next_chosen_node);

                travel_time = dist_matrix(current_node, next_chosen_node);
                arrival_at_chosen = departure_time_from_node + travel_time;
                service_start_chosen = max(arrival_at_chosen, all_nodes_early_tw(next_chosen_node));
                departure_time_from_node = service_start_chosen + all_nodes_service_time(next_chosen_node); % 更新离开当前服务节点的时间

                visited_customers(next_chosen_node-1) = true; % 标记该客户已被访问 (注意索引转换)
                current_node = next_chosen_node; % 更新当前节点
            end % 结束单条车辆路径构建

            % 路径末尾添加Depot，形成完整回路
            current_route = [current_route, 1];
            if length(current_route) > 2 % 只有当路径服务了至少一个客户时才添加
                ant_routes{end+1} = current_route;
            elseif ~all(visited_customers) && num_vehicles_ant > num_customers % 安全中断机制，防止死循环
                 warning('ACO: 尝试服务剩余客户时卡住。中断当前蚂蚁的构建过程。');
                 % 将所有剩余客户标记为已访问，以停止此蚂蚁
                 visited_customers(:) = true;
            end
        end % 结束 while ~all(visited_customers) (当前蚂蚁完成所有客户的分配)

        ant_solutions{ant} = ant_routes; % 存储当前蚂蚁的解
        if isempty(ant_routes) && num_customers > 0 % 如果没有形成任何路径但客户存在 (例如参数设置问题)
            ant_costs(ant) = inf; % 给予极大的成本惩罚
        else
            [cost_val, ~, ~] = calculate_total_cost(ant_routes, customer_data, depot_data, params);
            ant_costs(ant) = cost_val;
        end

    end % 结束蚂蚁循环 (所有蚂蚁均构建完成)

    % 找出当前迭代中的最优蚂蚁及其解
    [min_iter_cost, best_ant_idx] = min(ant_costs);
    iter_best_solution = ant_solutions{best_ant_idx};

    % 更新全局最优解
    if min_iter_cost < best_solution_cost
        best_solution_cost = min_iter_cost;
        best_solution_routes = iter_best_solution;
        fprintf('找到新的最优解: 成本 = %.2f\n', best_solution_cost);
    end
    best_costs_history(iter) = best_solution_cost; % 记录本次迭代的最优成本

    % 信息素更新 (蚁周模型 Ant Cycle Model)
    % 1. 信息素挥发
    tau = (1 - rho) * tau;

    % 2. 信息素增强 (所有蚂蚁根据其解的质量贡献信息素)
    for ant_idx = 1:num_ants
        if ~isinf(ant_costs(ant_idx)) && ant_costs(ant_idx) > 0 % 仅对有效解进行操作
            routes_this_ant = ant_solutions{ant_idx};
            delta_tau = Q / ant_costs(ant_idx); % 成本越低，增加的信息素越多
            for r = 1:length(routes_this_ant)
                route = routes_this_ant{r};
                for i = 1:(length(route)-1)
                    node1 = route(i);
                    node2 = route(i+1);
                    tau(node1, node2) = tau(node1, node2) + delta_tau;
                    tau(node2, node1) = tau(node2, node1) + delta_tau; % 对称路径信息素相同
                end
            end
        end
    end

    % (可选) 信息素限制 (Max-Min Ant System的特性, 非严格蚁周模型)
    % tau_min = 1e-4; tau_max = 10;
    % tau(tau < tau_min) = tau_min;
    % tau(tau > tau_max) = tau_max;

end % 结束迭代循环

fprintf('ACO算法结束。\n最优解成本: %.2f\n', best_solution_cost);
[final_cost, final_costs_breakdown, final_num_vehicles] = calculate_total_cost(best_solution_routes, customer_data, depot_data, params);
fprintf('最优解的成本构成:\n');
fprintf('  使用车辆数: %d\n', final_num_vehicles);
fprintf('  行驶成本: %.2f\n', final_costs_breakdown.distance);
fprintf('  车辆固定成本: %.2f\n', final_costs_breakdown.fixed_vehicle);
fprintf('  时间窗惩罚成本: %.2f\n', final_costs_breakdown.time_penalty);
fprintf('  货损成本: %.2f\n', final_costs_breakdown.cargo_damage);
fprintf('  温湿度控制成本: %.2f\n', final_costs_breakdown.refrigeration);
fprintf('  总成本 (验证): %.2f\n', final_cost);


%% 4. 可视化
% 迭代图: 显示算法收敛过程
figure;
plot(1:num_iterations, best_costs_history, 'LineWidth', 2);
xlabel('迭代次数');
ylabel('最优总成本');
title('CVRPTW的ACO算法收敛图');
grid on;

% 路径图: 显示最优解的车辆路径
figure;
hold on;
plot(depot_data.x, depot_data.y, 'ks', 'MarkerSize', 12, 'MarkerFaceColor', 'r'); % 绘制Depot
text(depot_data.x+1, depot_data.y+1, 'Depot');
plot([customer_data.x], [customer_data.y], 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); % 绘制客户点
for i = 1:num_customers
    text(customer_data(i).x+1, customer_data(i).y+1, sprintf('%d', i)); % 标记客户编号
end

colors = lines(length(best_solution_routes)); % 为不同路径生成不同颜色
if isempty(best_solution_routes)
    disp('未找到解决方案，无法绘制路径图。');
else
    for r = 1:length(best_solution_routes) % 遍历每条最优路径
        route = best_solution_routes{r};
        route_coords_x = all_nodes_coords(route, 1); % 获取路径中各节点的x坐标
        route_coords_y = all_nodes_coords(route, 2); % 获取路径中各节点的y坐标
        plot(route_coords_x, route_coords_y, 'LineWidth', 1.5, 'Color', colors(r,:), 'Marker', '>'); % 绘制路径
    end
end

axis equal; % 保持x,y轴比例一致
title(['找到的最优路径 (成本: ', sprintf('%.2f', best_solution_cost), ')']);
xlabel('X 坐标');
ylabel('Y 坐标');
legend('Depot', '客户点', 'Location', 'bestoutside'); % 图例
hold off;

% 在路径规划完成后，输出每辆车的路径
if isempty(solution_routes) || isinf(solution_cost)
    fprintf('未找到有效路径。\n');
else
    fprintf('所有车辆的路径规划如下：\n');
    for r = 1:length(solution_routes)
        route = solution_routes{r};
        fprintf('车辆 %d 的路径: ', r);
        % 遍历并打印路径上的节点编号
        for i = 1:length(route)
            if i == length(route)
                fprintf('%d\n', route(i)); % 最后一个节点换行
            else
                fprintf('%d -> ', route(i)); % 连接符号
            end
        end
    end
end

rng(s); % 恢复随机数生成器的状态