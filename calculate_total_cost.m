function [total_cost, costs_breakdown, num_vehicles] = calculate_total_cost(routes, customer_data, depot_data, params)
    % 计算给定解决方案的总成本
    % 输入:
    %   routes: 路径的单元格数组, 例如: {[1 2 3 1], [1 4 5 1]} (1代表depot)
    %   customer_data: 客户数据的结构体数组，包含字段: x, y, demand, early_tw, late_tw, service_time, cargo_value, damage_rate
    %   depot_data: 配送中心数据的结构体，包含字段: x, y, early_tw, late_tw
    %   params: 成本参数的结构体: vehicle_capacity, fixed_vehicle_cost, cost_per_distance,
    %                               penalty_per_late_time, unit_refrigeration_cost
    % 输出:
    %   total_cost: 总成本
    %   costs_breakdown: 各项成本明细的结构体
    %   num_vehicles: 使用的车辆数量

    num_nodes = length(customer_data) + 1; % 总节点数 (包括depot作为节点1)
    % 将depot和客户数据整合，方便索引
    % 节点1是depot, 节点2到N+1是客户
    all_nodes_coords = [depot_data.x, depot_data.y; [customer_data.x]', [customer_data.y]']; % 所有节点坐标
    all_nodes_demand = [0; [customer_data.demand]']; % 所有节点需求量 (depot需求为0)
    all_nodes_early_tw = [depot_data.early_tw; [customer_data.early_tw]']; % 所有节点最早到达时间
    all_nodes_late_tw = [depot_data.late_tw; [customer_data.late_tw]'];   % 所有节点最晚到达时间
    all_nodes_service_time = [0; [customer_data.service_time]']; % 所有节点服务时长 (depot服务时长为0)

    total_distance_cost = 0;        % 总行驶成本
    total_fixed_vehicle_cost = 0;   % 总车辆固定成本
    total_time_penalty_cost = 0;    % 总时间窗惩罚成本
    total_cargo_damage_cost = 0;    % 总货损成本
    total_refrigeration_cost = 0;   % 总制冷成本

    num_vehicles = length(routes); % 使用的车辆数等于路径数
    total_fixed_vehicle_cost = num_vehicles * params.fixed_vehicle_cost;

    for r = 1:num_vehicles % 遍历每条路径 (即每辆车)
        route = routes{r};
        current_time = all_nodes_early_tw(1); % 车辆从depot出发的初始时间，默认为depot的最早可出发时间
        current_load = 0;                     % 当前车辆载重
        route_distance = 0;                   % 当前路径行驶距离
        % route_refrigeration_duration = 0;  % 当前路径制冷时长 (在下方精确计算)

        % 车辆从上一个节点出发的时间
        departure_time_from_prev_node = current_time;

        for k = 1:(length(route) - 1) % 遍历路径中的每一段 (node_from -> node_to)
            node_from = route(k);
            node_to = route(k+1);

            dist_segment = norm(all_nodes_coords(node_from,:) - all_nodes_coords(node_to,:)); % 段距离
            route_distance = route_distance + dist_segment;
            travel_time = dist_segment; % 假设车速为1 单位距离/单位时间

            arrival_time_at_node_to = departure_time_from_prev_node + travel_time; % 到达 node_to 的时间

            if node_to ~= 1 % 如果不是返回depot
                % 等待时间计算 (早到需等待)
                wait_time = max(0, all_nodes_early_tw(node_to) - arrival_time_at_node_to);
                service_start_time = arrival_time_at_node_to + wait_time; % 服务开始时间

                % 晚到惩罚计算 (仅晚于最晚服务时间才惩罚)
                lateness = max(0, service_start_time - all_nodes_late_tw(node_to));
                total_time_penalty_cost = total_time_penalty_cost + lateness * params.penalty_per_late_time;

                % 计算该客户的货损成本
                % customer_data 的索引是 node_to - 1 (因为customer_data是从1开始的客户索引，而node_to中depot是1)
                customer_idx_in_data = node_to - 1;
                total_cargo_damage_cost = total_cargo_damage_cost + ...
                                          customer_data(customer_idx_in_data).demand * ...
                                          customer_data(customer_idx_in_data).cargo_value * ...
                                          customer_data(customer_idx_in_data).damage_rate;

                % 更新离开当前节点的时间 (服务完成后)
                departure_time_from_prev_node = service_start_time + all_nodes_service_time(node_to);
                current_load = current_load + all_nodes_demand(node_to); % 更新载重 (这应该在路径构建时检查容量)
            else % 如果是返回depot
                departure_time_from_prev_node = arrival_time_at_node_to; % 返回depot没有服务时间
                % 检查返回depot是否在其时间窗内 (如果有严格约束)
                 lateness_depot = max(0, arrival_time_at_node_to - all_nodes_late_tw(1)); % depot也可能有最晚返回时间
                 total_time_penalty_cost = total_time_penalty_cost + lateness_depot * params.penalty_per_late_time;
            end
            current_time = departure_time_from_prev_node; % 更新当前时间为此节点出发时间，用于下一段路程
        end

        % 计算当前路径的温湿度控制成本
        % 制冷时间等于车辆离开depot到返回depot的总时长
        if ~isempty(route) && route(1) == 1 && route(end) == 1 && length(route) > 1
            % 为了精确计算路径的总运行时长（设备时间），需要重新追溯时间
            route_start_time_for_refrigeration = all_nodes_early_tw(1); % 假设每辆车从depot的最早时间出发
            time_ptr = route_start_time_for_refrigeration;
            prev_node_in_trace = route(1); % 起始是depot

            for k_trace = 2:length(route) % 从路径的第二个节点开始（第一个客户或返回的depot）
                curr_node_in_trace = route(k_trace);
                dist_seg_trace = norm(all_nodes_coords(prev_node_in_trace,:) - all_nodes_coords(curr_node_in_trace,:));
                travel_t_trace = dist_seg_trace; % 假设速度为1

                arrival_t_trace = time_ptr + travel_t_trace; % 到达当前追踪节点的时间

                if curr_node_in_trace ~= 1 % 如果不是depot
                    service_start_t_trace = max(arrival_t_trace, all_nodes_early_tw(curr_node_in_trace)); % 服务开始时间
                    time_ptr = service_start_t_trace + all_nodes_service_time(curr_node_in_trace);      % 服务结束离开的时间
                else % 如果是返回depot
                    time_ptr = arrival_t_trace; % 到达depot的时间
                end
                prev_node_in_trace = curr_node_in_trace;
            end
            route_refrigeration_duration = time_ptr - route_start_time_for_refrigeration; % 总制冷时长
        else
            route_refrigeration_duration = 0; % 对于无效路径或空路径
        end

        total_refrigeration_cost = total_refrigeration_cost + route_refrigeration_duration * params.unit_refrigeration_cost;
        total_distance_cost = total_distance_cost + route_distance * params.cost_per_distance;
    end

    % 计算总成本
    total_cost = total_distance_cost + total_fixed_vehicle_cost + ...
                 total_time_penalty_cost + total_cargo_damage_cost + total_refrigeration_cost;

    % 存储各项成本明细
    costs_breakdown.distance = total_distance_cost;
    costs_breakdown.fixed_vehicle = total_fixed_vehicle_cost;
    costs_breakdown.time_penalty = total_time_penalty_cost;
    costs_breakdown.cargo_damage = total_cargo_damage_cost;
    costs_breakdown.refrigeration = total_refrigeration_cost;
end