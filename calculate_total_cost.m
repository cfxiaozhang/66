function [total_cost, costs_breakdown, num_vehicles] = calculate_total_cost(routes, customer_data, depot_data, params)
    % �����������������ܳɱ�
    % ����:
    %   routes: ·���ĵ�Ԫ������, ����: {[1 2 3 1], [1 4 5 1]} (1����depot)
    %   customer_data: �ͻ����ݵĽṹ�����飬�����ֶ�: x, y, demand, early_tw, late_tw, service_time, cargo_value, damage_rate
    %   depot_data: �����������ݵĽṹ�壬�����ֶ�: x, y, early_tw, late_tw
    %   params: �ɱ������Ľṹ��: vehicle_capacity, fixed_vehicle_cost, cost_per_distance,
    %                               penalty_per_late_time, unit_refrigeration_cost
    % ���:
    %   total_cost: �ܳɱ�
    %   costs_breakdown: ����ɱ���ϸ�Ľṹ��
    %   num_vehicles: ʹ�õĳ�������

    num_nodes = length(customer_data) + 1; % �ܽڵ��� (����depot��Ϊ�ڵ�1)
    % ��depot�Ϳͻ��������ϣ���������
    % �ڵ�1��depot, �ڵ�2��N+1�ǿͻ�
    all_nodes_coords = [depot_data.x, depot_data.y; [customer_data.x]', [customer_data.y]']; % ���нڵ�����
    all_nodes_demand = [0; [customer_data.demand]']; % ���нڵ������� (depot����Ϊ0)
    all_nodes_early_tw = [depot_data.early_tw; [customer_data.early_tw]']; % ���нڵ����絽��ʱ��
    all_nodes_late_tw = [depot_data.late_tw; [customer_data.late_tw]'];   % ���нڵ�������ʱ��
    all_nodes_service_time = [0; [customer_data.service_time]']; % ���нڵ����ʱ�� (depot����ʱ��Ϊ0)

    total_distance_cost = 0;        % ����ʻ�ɱ�
    total_fixed_vehicle_cost = 0;   % �ܳ����̶��ɱ�
    total_time_penalty_cost = 0;    % ��ʱ�䴰�ͷ��ɱ�
    total_cargo_damage_cost = 0;    % �ܻ���ɱ�
    total_refrigeration_cost = 0;   % ������ɱ�

    num_vehicles = length(routes); % ʹ�õĳ���������·����
    total_fixed_vehicle_cost = num_vehicles * params.fixed_vehicle_cost;

    for r = 1:num_vehicles % ����ÿ��·�� (��ÿ����)
        route = routes{r};
        current_time = all_nodes_early_tw(1); % ������depot�����ĳ�ʼʱ�䣬Ĭ��Ϊdepot������ɳ���ʱ��
        current_load = 0;                     % ��ǰ��������
        route_distance = 0;                   % ��ǰ·����ʻ����
        % route_refrigeration_duration = 0;  % ��ǰ·������ʱ�� (���·���ȷ����)

        % ��������һ���ڵ������ʱ��
        departure_time_from_prev_node = current_time;

        for k = 1:(length(route) - 1) % ����·���е�ÿһ�� (node_from -> node_to)
            node_from = route(k);
            node_to = route(k+1);

            dist_segment = norm(all_nodes_coords(node_from,:) - all_nodes_coords(node_to,:)); % �ξ���
            route_distance = route_distance + dist_segment;
            travel_time = dist_segment; % ���賵��Ϊ1 ��λ����/��λʱ��

            arrival_time_at_node_to = departure_time_from_prev_node + travel_time; % ���� node_to ��ʱ��

            if node_to ~= 1 % ������Ƿ���depot
                % �ȴ�ʱ����� (�絽��ȴ�)
                wait_time = max(0, all_nodes_early_tw(node_to) - arrival_time_at_node_to);
                service_start_time = arrival_time_at_node_to + wait_time; % ����ʼʱ��

                % ���ͷ����� (�������������ʱ��ųͷ�)
                lateness = max(0, service_start_time - all_nodes_late_tw(node_to));
                total_time_penalty_cost = total_time_penalty_cost + lateness * params.penalty_per_late_time;

                % ����ÿͻ��Ļ���ɱ�
                % customer_data �������� node_to - 1 (��Ϊcustomer_data�Ǵ�1��ʼ�Ŀͻ���������node_to��depot��1)
                customer_idx_in_data = node_to - 1;
                total_cargo_damage_cost = total_cargo_damage_cost + ...
                                          customer_data(customer_idx_in_data).demand * ...
                                          customer_data(customer_idx_in_data).cargo_value * ...
                                          customer_data(customer_idx_in_data).damage_rate;

                % �����뿪��ǰ�ڵ��ʱ�� (������ɺ�)
                departure_time_from_prev_node = service_start_time + all_nodes_service_time(node_to);
                current_load = current_load + all_nodes_demand(node_to); % �������� (��Ӧ����·������ʱ�������)
            else % ����Ƿ���depot
                departure_time_from_prev_node = arrival_time_at_node_to; % ����depotû�з���ʱ��
                % ��鷵��depot�Ƿ�����ʱ�䴰�� (������ϸ�Լ��)
                 lateness_depot = max(0, arrival_time_at_node_to - all_nodes_late_tw(1)); % depotҲ������������ʱ��
                 total_time_penalty_cost = total_time_penalty_cost + lateness_depot * params.penalty_per_late_time;
            end
            current_time = departure_time_from_prev_node; % ���µ�ǰʱ��Ϊ�˽ڵ����ʱ�䣬������һ��·��
        end

        % ���㵱ǰ·������ʪ�ȿ��Ƴɱ�
        % ����ʱ����ڳ����뿪depot������depot����ʱ��
        if ~isempty(route) && route(1) == 1 && route(end) == 1 && length(route) > 1
            % Ϊ�˾�ȷ����·����������ʱ�����豸ʱ�䣩����Ҫ����׷��ʱ��
            route_start_time_for_refrigeration = all_nodes_early_tw(1); % ����ÿ������depot������ʱ�����
            time_ptr = route_start_time_for_refrigeration;
            prev_node_in_trace = route(1); % ��ʼ��depot

            for k_trace = 2:length(route) % ��·���ĵڶ����ڵ㿪ʼ����һ���ͻ��򷵻ص�depot��
                curr_node_in_trace = route(k_trace);
                dist_seg_trace = norm(all_nodes_coords(prev_node_in_trace,:) - all_nodes_coords(curr_node_in_trace,:));
                travel_t_trace = dist_seg_trace; % �����ٶ�Ϊ1

                arrival_t_trace = time_ptr + travel_t_trace; % ���ﵱǰ׷�ٽڵ��ʱ��

                if curr_node_in_trace ~= 1 % �������depot
                    service_start_t_trace = max(arrival_t_trace, all_nodes_early_tw(curr_node_in_trace)); % ����ʼʱ��
                    time_ptr = service_start_t_trace + all_nodes_service_time(curr_node_in_trace);      % ��������뿪��ʱ��
                else % ����Ƿ���depot
                    time_ptr = arrival_t_trace; % ����depot��ʱ��
                end
                prev_node_in_trace = curr_node_in_trace;
            end
            route_refrigeration_duration = time_ptr - route_start_time_for_refrigeration; % ������ʱ��
        else
            route_refrigeration_duration = 0; % ������Ч·�����·��
        end

        total_refrigeration_cost = total_refrigeration_cost + route_refrigeration_duration * params.unit_refrigeration_cost;
        total_distance_cost = total_distance_cost + route_distance * params.cost_per_distance;
    end

    % �����ܳɱ�
    total_cost = total_distance_cost + total_fixed_vehicle_cost + ...
                 total_time_penalty_cost + total_cargo_damage_cost + total_refrigeration_cost;

    % �洢����ɱ���ϸ
    costs_breakdown.distance = total_distance_cost;
    costs_breakdown.fixed_vehicle = total_fixed_vehicle_cost;
    costs_breakdown.time_penalty = total_time_penalty_cost;
    costs_breakdown.cargo_damage = total_cargo_damage_cost;
    costs_breakdown.refrigeration = total_refrigeration_cost;
end