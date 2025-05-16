clear; clc; close all; % 清理工作区，命令窗口，关闭所有图形

%% 1. 问题定义 (与ACO版本相同)
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

for i=1:num_customers
    customer_data(i).late_tw = max(customer_data(i).late_tw, customer_data(i).early_tw + customer_data(i).service_time + 10);
end

depot_data.x = 35;
depot_data.y = 35;
depot_data.early_tw = 0;
depot_data.late_tw = 230; % 仓库最晚关闭时间

params.vehicle_capacity = 100;
params.fixed_vehicle_cost = 300;
params.cost_per_distance = 3;
params.penalty_per_late_time = 3;
params.unit_refrigeration_cost = 1;

all_nodes_coords = [depot_data.x, depot_data.y; [customer_data.x]', [customer_data.y]'];
all_nodes_demand = [0; [customer_data.demand]'];
all_nodes_early_tw = [depot_data.early_tw; [customer_data.early_tw]'];
all_nodes_late_tw = [depot_data.late_tw; [customer_data.late_tw]'];
all_nodes_service_time = [0; [customer_data.service_time]'];

num_total_nodes = size(all_nodes_coords, 1);

dist_matrix = zeros(num_total_nodes, num_total_nodes);
for i = 1:num_total_nodes
    for j = i+1:num_total_nodes
        dist_matrix(i,j) = norm(all_nodes_coords(i,:) - all_nodes_coords(j,:));
        dist_matrix(j,i) = dist_matrix(i,j);
    end
end

%% 2. 最近邻算法 (NN) 实现
fprintf('CVRPTW的最近邻算法开始...\n');

solution_routes = {};       % 存储最终解的路径
visited_customers = false(1, num_customers); % 标记客户是否已被访问 (客户索引1到num_customers)
num_vehicles_used = 0;

while ~all(visited_customers) % 当仍有客户未被服务时
    num_vehicles_used = num_vehicles_used + 1;
    current_route = [1]; % 新路径从Depot (节点1) 开始
    current_load = 0;    % 当前车辆载重
    current_node = 1;    % 当前所在节点
    % 车辆从仓库出发的初始时间，等于仓库最早可出发时间
    departure_time_from_current_node = all_nodes_early_tw(1); 

    while true % 构建单条车辆路径
        best_next_node = -1;          % 存储找到的最佳下一个节点
        min_distance_to_next = inf;   % 到最佳下一个节点的最小距离
        
        % 为最佳下一个节点预计算的时间信息 (如果被选中)
        arrival_at_best_next_cand = -1;
        service_start_at_best_next_cand = -1;
        departure_from_best_next_cand = -1;

        % 寻找下一个可服务的、最近的客户节点
        for candidate_node_idx = 2:num_total_nodes % 遍历所有客户节点 (索引2到N+1)
            customer_actual_idx = candidate_node_idx - 1; % 对应 customer_data 和 visited_customers 中的索引

            if ~visited_customers(customer_actual_idx) && ... % 如果该客户未被全局访问过
               (current_load + all_nodes_demand(candidate_node_idx) <= params.vehicle_capacity) % 且不超过车辆容量

                % 时间可行性检查
                travel_time_to_candidate = dist_matrix(current_node, candidate_node_idx);
                arrival_at_candidate = departure_time_from_current_node + travel_time_to_candidate;
                
                % 早到需等待，计算服务开始时间
                service_start_at_candidate = max(arrival_at_candidate, all_nodes_early_tw(candidate_node_idx));
                
                % 离开候选客户点的时间
                departure_from_candidate = service_start_at_candidate + all_nodes_service_time(candidate_node_idx);
                
                % 从候选客户点返回仓库所需时间
                time_to_return_depot_from_candidate = dist_matrix(candidate_node_idx, 1);
                
                % 检查：服务完候选客户并返回仓库，是否会超过仓库的最晚关闭时间
                if departure_from_candidate + time_to_return_depot_from_candidate <= all_nodes_late_tw(1) % Depot最晚时间
                    % 如果满足容量和时间约束，则考虑距离
                    distance_to_candidate = dist_matrix(current_node, candidate_node_idx);
                    if distance_to_candidate < min_distance_to_next
                        min_distance_to_next = distance_to_candidate;
                        best_next_node = candidate_node_idx;
                        % 记录这个候选者的时间信息，如果它最终被选为最佳
                        arrival_at_best_next_cand = arrival_at_candidate;
                        service_start_at_best_next_cand = service_start_at_candidate;
                        departure_from_best_next_cand = departure_from_candidate;
                    end
                end
            end
        end % 结束对候选节点的遍历

        if best_next_node == -1 % 如果没有找到可行的下一个客户
            break; % 结束当前车辆的路径构建
        else
            % 将找到的最佳下一个节点添加到当前路径
            current_route = [current_route, best_next_node];
            current_load = current_load + all_nodes_demand(best_next_node);
            visited_customers(best_next_node-1) = true; % 标记该客户已被访问 (注意索引转换)
            current_node = best_next_node; % 更新当前节点
            
            % 更新从当前节点（即新加入的best_next_node）出发的时间
            % 这个时间是服务完best_next_node后的离开时间
            departure_time_from_current_node = departure_from_best_next_cand;
        end
    end % 结束单条车辆路径构建

    % 路径末尾添加Depot，形成完整回路
    current_route = [current_route, 1];
    if length(current_route) > 2 % 只有当路径服务了至少一个客户时才是一条有效路径
        solution_routes{end+1} = current_route;
    elseif ~all(visited_customers) && num_vehicles_used > num_customers % 安全中断，防止所有客户都无法服务时死循环
        fprintf('警告: 无法为所有剩余客户规划路径，部分客户可能未被服务。\n');
        break; % 退出主循环
    end
    
    if isempty(solution_routes) && ~all(visited_customers) && num_customers > 0
         % 如果一条路径都建不起来，但还有客户未服务，说明可能有问题
         % 比如第一个客户就无法满足约束
         fprintf('警告: 未能构建任何有效路径，但仍有客户未服务。\n');
         break;
    end

end % 结束 while ~all(visited_customers)

%% 3. 计算总成本并输出结果
if ~all(visited_customers) && num_customers > 0
    fprintf('注意: 并非所有客户都被服务!\n');
    fprintf('已服务客户数: %d/%d\n', sum(visited_customers), num_customers);
end

if isempty(solution_routes) && num_customers > 0
    fprintf('最近邻算法未能找到任何有效路径。\n');
    solution_cost = inf;
    costs_breakdown = struct('distance', inf, 'fixed_vehicle', inf, 'time_penalty', inf, 'cargo_damage', inf, 'refrigeration', inf);
    actual_num_vehicles = 0;
else
    [solution_cost, costs_breakdown, actual_num_vehicles] = calculate_total_cost(solution_routes, customer_data, depot_data, params);
end

fprintf('最近邻算法结束。\n');
fprintf('找到的解的成本: %.2f\n', solution_cost);
if ~isinf(solution_cost)
    fprintf('成本构成:\n');
    fprintf('  使用车辆数: %d\n', actual_num_vehicles);
    fprintf('  行驶成本: %.2f\n', costs_breakdown.distance);
    fprintf('  车辆固定成本: %.2f\n', costs_breakdown.fixed_vehicle);
    fprintf('  时间窗惩罚成本: %.2f\n', costs_breakdown.time_penalty);
    fprintf('  货损成本: %.2f\n', costs_breakdown.cargo_damage);
    fprintf('  温湿度控制成本: %.2f\n', costs_breakdown.refrigeration);
end

%% 4. 可视化
% "迭代图" 对于NN意义不大，因为它只产生一个解。这里画一个点表示最终成本。
figure;
plot(1, solution_cost, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlim([0 2]); % 调整X轴范围以便显示单个点
xlabel('解 (最近邻算法)');
ylabel('总成本');
title('最近邻算法求解CVRPTW的成本');
grid on;
if isinf(solution_cost)
    text(1, 0, '无解或成本无限大', 'HorizontalAlignment', 'center');
else
    text(1, solution_cost, sprintf('%.2f', solution_cost), 'VerticalAlignment','bottom', 'HorizontalAlignment','center');
end


% 路径图
figure;
hold on;
plot(depot_data.x, depot_data.y, 'ks', 'MarkerSize', 12, 'MarkerFaceColor', 'r'); % 绘制Depot
text(depot_data.x+1, depot_data.y+1, 'Depot');
plot([customer_data.x], [customer_data.y], 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); % 绘制客户点
for i = 1:num_customers
    text(customer_data(i).x+1, customer_data(i).y+1, sprintf('%d', i)); % 标记客户编号
end

if isempty(solution_routes) || isinf(solution_cost)
    disp('无有效路径可供绘制。');
else
    colors = lines(length(solution_routes)); % 为不同路径生成不同颜色
    for r = 1:length(solution_routes) % 遍历每条路径
        route = solution_routes{r};
        route_coords_x = all_nodes_coords(route, 1); % 获取路径中各节点的x坐标
        route_coords_y = all_nodes_coords(route, 2); % 获取路径中各节点的y坐标
        plot(route_coords_x, route_coords_y, 'LineWidth', 1.5, 'Color', colors(r,:), 'Marker', '>'); % 绘制路径
    end
end

axis equal; % 保持x,y轴比例一致
if isinf(solution_cost)
    title('最近邻算法路径图 (无有效解)');
else
    title(['最近邻算法找到的路径 (成本: ', sprintf('%.2f', solution_cost), ')']);
end
xlabel('X 坐标');
ylabel('Y 坐标');
legend('Depot', '客户点', 'Location', 'bestoutside'); % 图例
hold off;

rng(s); % 恢复随机数生成器的状态