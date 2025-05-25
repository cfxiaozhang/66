clear; clc; close all; % ��������������ڣ��ر�����ͼ��

%% 1. ���ⶨ��
% �ͻ�����: x����, y����, ������, ���絽��ʱ��, ������ʱ��, ����ʱ��, ���ﵥλ��ֵ, ������
s = rng(1); % ��������������Ա����ɸ���

num_customers = 50; % �ͻ�����
customer_data = struct('x', num2cell(xlsread('�ͻ���.xlsx','A2:AX2')), ... % x����
                       'y', num2cell(xlsread('�ͻ���.xlsx','A3:AX3')), ... % y����
                       'demand', num2cell(xlsread('�ͻ���.xlsx','A4:AX4')), ... % ������
                       'early_tw', num2cell(xlsread('�ͻ���.xlsx','A5:AX5')), ... % ���絽��ʱ��
                       'late_tw', num2cell(xlsread('�ͻ���.xlsx','A6:AX6')), ...% ������ʱ��
                       'service_time', num2cell(xlsread('�ͻ���.xlsx','A7:AX7')), ... % ����ʱ��
                       'cargo_value', num2cell(1000 * ones(1, num_customers)), ...    % ���ﵥλ��ֵ
                       'damage_rate', num2cell(0.001 * ones(1, num_customers)));   % ������

% ȷ������ʱ�䴰���� (����ʱ�䴰 + ����ʱ��)��������ʱ�䴰�ᴦ�����������ø�����
for i=1:num_customers
    customer_data(i).late_tw = max(customer_data(i).late_tw, customer_data(i).early_tw + customer_data(i).service_time + 10); % ��10��Ϊ����Щ���
end

% �������� (Depot) ����
depot_data.x = 35;         % Depot x����
depot_data.y = 35;         % Depot y����
depot_data.early_tw = 0;   % Depot ����ɳ���ʱ��
depot_data.late_tw = 230;  % Depot ������Ӫ����ʱ�� (���������ڴ�֮ǰ����)

% ��������
params.vehicle_capacity = 100;         % ��������
params.fixed_vehicle_cost = 300;      % ÿ�����Ĺ̶��ɱ�
params.cost_per_distance = 3;         % ��λ�������ʻ�ɱ�
params.penalty_per_late_time = 3;     % ��λ��ʱ��ĳͷ��ɱ�
params.unit_refrigeration_cost = 1; % ��λ��ʪ�ȿ���ʱ��ĳɱ�

% ����Depot�Ϳͻ����ݣ�������������������
% �ڵ�1: Depot, �ڵ� 2 �� N+1: �ͻ�
all_nodes_coords = [depot_data.x, depot_data.y; [customer_data.x]', [customer_data.y]']; % ���нڵ�����
all_nodes_demand = [0; [customer_data.demand]']; % ���нڵ����� (DepotΪ0)
all_nodes_early_tw = [depot_data.early_tw; [customer_data.early_tw]']; % ���нڵ�����ʱ�䴰
all_nodes_late_tw = [depot_data.late_tw; [customer_data.late_tw]'];   % ���нڵ�����ʱ�䴰
all_nodes_service_time = [0; [customer_data.service_time]']; % ���нڵ����ʱ�� (DepotΪ0)

num_total_nodes = size(all_nodes_coords, 1); % �ܽڵ���

% ����������
dist_matrix = zeros(num_total_nodes, num_total_nodes);
for i = 1:num_total_nodes
    for j = i+1:num_total_nodes
        dist_matrix(i,j) = norm(all_nodes_coords(i,:) - all_nodes_coords(j,:));
        dist_matrix(j,i) = dist_matrix(i,j);
    end
end

%% 2. ��Ⱥ�㷨 (ACO) ����
num_ants = 100;        % ��������
num_iterations = 300; % ��������
alpha = 1;            % ��Ϣ����Ҫ�̶�����
beta = 3;             % ����ʽ��Ϣ (����) ��Ҫ�̶�����
rho = 0.5;            % ��Ϣ�ػӷ���
Q = 100;              % ��Ϣ������ǿ��ϵ��

% ����ʽ��Ϣ (eta = 1/����)
eta = 1 ./ dist_matrix;
eta(isinf(eta) | isnan(eta)) = 0; % ������������� (ͬһ���ڵ����Ϊ0)

% ��ʼ����Ϣ�ؾ��� (tau)
tau = ones(num_total_nodes, num_total_nodes) * 1e-3; % ��ʼ��Ϣ����Ϊһ����С��ֵ

best_solution_routes = {};         % �洢ȫ�����Ž��·��
best_solution_cost = inf;          % �洢ȫ�����Ž�ĳɱ�
best_costs_history = zeros(1, num_iterations); % ��¼ÿ�ε��������ųɱ�

%% 3. ACO ��ѭ��
fprintf('VRPTW����Ⱥ�㷨��ʼ...\n');
for iter = 1:num_iterations
    fprintf('���� %d/%d\n', iter, num_iterations);

    ant_solutions = cell(num_ants, 1); % �洢��ǰ������ÿֻ���ϵĽ� (·������)
    ant_costs = zeros(num_ants, 1);    % �洢��ǰ������ÿֻ���Ͻ�ĳɱ�

    for ant = 1:num_ants % ����ÿֻ����
        ant_routes = {}; % ��ǰ���Ϲ�����·������
        visited_customers = false(1, num_customers); % ��ǿͻ��Ƿ��ѱ����� (�ͻ�������1��num_customers, ��Ӧ�ڵ�����2��N+1)
        num_vehicles_ant = 0; % ��ǰ����ʹ�õĳ�����

        while ~all(visited_customers) % �����пͻ�δ������ʱ
            num_vehicles_ant = num_vehicles_ant + 1;
            current_route = [1]; % ��·����Depot (�ڵ�1) ��ʼ
            current_load = 0;    % ��ǰ��������
            current_time = all_nodes_early_tw(1); % ������Depot������ʱ��
            current_node = 1;    % ��ǰ���ڽڵ�

            departure_time_from_node = current_time; % �ӵ�ǰ�ڵ������ʱ��

            while true % ������������·��
                possible_next_nodes = []; % ���ܵ���һ���ͻ��ڵ��б�
                probabilities = [];       % ��Ӧ�ڵ��ѡ�����

                % Ѱ����һ���ɷ���Ŀͻ��ڵ�
                for next_node_idx = 2:num_total_nodes % �������пͻ��ڵ� (����2��N+1)
                    customer_actual_idx = next_node_idx - 1; % ��Ӧ customer_data �е�����
                    if ~visited_customers(customer_actual_idx) && ... % ����ÿͻ�δ��ȫ�ַ��ʹ�
                       (current_load + all_nodes_demand(next_node_idx) <= params.vehicle_capacity) % �Ҳ�������������

                        % ʱ�䴰�����Գ������ (����ϸ�ļ���ڳɱ����㺯����)
                        travel_time_to_next = dist_matrix(current_node, next_node_idx);
                        arrival_at_next = departure_time_from_node + travel_time_to_next;

                        % ����һ���ڵ���������ʼʱ��
                        service_start_at_next = max(arrival_at_next, all_nodes_early_tw(next_node_idx));

                        % ������һ���ڵ������Ϻ󷵻�depot�Ƿ�ᳬ��depot��������Ӫʱ��
                        departure_from_next = service_start_at_next + all_nodes_service_time(next_node_idx);
                        time_to_return_depot = dist_matrix(next_node_idx, 1); % ��next_node_idx����depot��ʱ��

                        if departure_from_next + time_to_return_depot <= all_nodes_late_tw(1) % Depot������ر�ʱ��
                             possible_next_nodes = [possible_next_nodes, next_node_idx];
                        end
                    end
                end

                if isempty(possible_next_nodes) % ���û�п�ѡ�����һ���ڵ�
                    break; % ������ǰ������·������ (���ܳ�����������ʣ��ͻ�������ʱ��Լ��)
                end

                % ����ѡ����һ���ڵ�ĸ���
                prob_numerator = []; % �洢���ʹ�ʽ�ķ��Ӳ���
                for nn_idx = 1:length(possible_next_nodes)
                    node_cand = possible_next_nodes(nn_idx); % ��ѡ�ڵ�
                    prob_numerator(nn_idx) = (tau(current_node, node_cand)^alpha) * ... % ��Ϣ��Ũ��
                                             (eta(current_node, node_cand)^beta);      % ����ʽ��Ϣ
                end

                if sum(prob_numerator) == 0 % ��������� (������з��Ӷ�Ϊ0)
                    % ����������ܷ�����etaΪ0 (����beta�ܸ��Ҿ���ܴ�) ��tau��С��
                    % �򵥴������ѡ��һ��������ֹ��ǰ����·����
                     if isempty(possible_next_nodes)
                        break;
                     else % ����ж������Ϊ0��ѡ����ѡһ��
                        selected_idx_in_possible = randi(length(possible_next_nodes));
                     end
                else
                    probabilities = prob_numerator / sum(prob_numerator); % ��һ���õ�����
                    % ���̶�ѡ����һ���ڵ�
                    cumulative_prob = cumsum(probabilities);
                    r_select = rand();
                    selected_idx_in_possible = find(r_select <= cumulative_prob, 1, 'first');
                end

                if isempty(selected_idx_in_possible) % ���û��ѡ���ڵ� (�����ϲ�Ӧ����������possible_next_nodesΪ�ջ�sum(prob_num)<=0)
                     break; % ��ֹ��ǰ����·��
                end

                next_chosen_node = possible_next_nodes(selected_idx_in_possible); % ѡ�е���һ���ڵ�

                % ����·�������ء�ʱ�����Ϣ
                current_route = [current_route, next_chosen_node];
                current_load = current_load + all_nodes_demand(next_chosen_node);

                travel_time = dist_matrix(current_node, next_chosen_node);
                arrival_at_chosen = departure_time_from_node + travel_time;
                service_start_chosen = max(arrival_at_chosen, all_nodes_early_tw(next_chosen_node));
                departure_time_from_node = service_start_chosen + all_nodes_service_time(next_chosen_node); % �����뿪��ǰ����ڵ��ʱ��

                visited_customers(next_chosen_node-1) = true; % ��Ǹÿͻ��ѱ����� (ע������ת��)
                current_node = next_chosen_node; % ���µ�ǰ�ڵ�
            end % ������������·������

            % ·��ĩβ���Depot���γ�������·
            current_route = [current_route, 1];
            if length(current_route) > 2 % ֻ�е�·������������һ���ͻ�ʱ�����
                ant_routes{end+1} = current_route;
            elseif ~all(visited_customers) && num_vehicles_ant > num_customers % ��ȫ�жϻ��ƣ���ֹ��ѭ��
                 warning('ACO: ���Է���ʣ��ͻ�ʱ��ס���жϵ�ǰ���ϵĹ������̡�');
                 % ������ʣ��ͻ����Ϊ�ѷ��ʣ���ֹͣ������
                 visited_customers(:) = true;
            end
        end % ���� while ~all(visited_customers) (��ǰ����������пͻ��ķ���)

        ant_solutions{ant} = ant_routes; % �洢��ǰ���ϵĽ�
        if isempty(ant_routes) && num_customers > 0 % ���û���γ��κ�·�����ͻ����� (���������������)
            ant_costs(ant) = inf; % ���輫��ĳɱ��ͷ�
        else
            [cost_val, ~, ~] = calculate_total_cost(ant_routes, customer_data, depot_data, params);
            ant_costs(ant) = cost_val;
        end

    end % ��������ѭ�� (�������Ͼ��������)

    % �ҳ���ǰ�����е��������ϼ����
    [min_iter_cost, best_ant_idx] = min(ant_costs);
    iter_best_solution = ant_solutions{best_ant_idx};

    % ����ȫ�����Ž�
    if min_iter_cost < best_solution_cost
        best_solution_cost = min_iter_cost;
        best_solution_routes = iter_best_solution;

        % ������ϸ�ɱ�
        [new_cost, new_costs_breakdown, new_num_vehicles] = calculate_total_cost(best_solution_routes, customer_data, depot_data, params);

        % ��ӡ�µ����Ž����ɱ�
        fprintf('\n�ҵ��µ����Ž� (���� %d):\n', iter);
        fprintf('  �ܳɱ�: %.2f\n', new_cost);
        fprintf('  ʹ�ó�����: %d\n', new_num_vehicles);
        fprintf('  ����ɱ�: %.2f\n', new_costs_breakdown.distance);
        fprintf('  �̶������ɱ�: %.2f\n', new_costs_breakdown.fixed_vehicle);
        fprintf('  ʱ�䴰����: %.2f\n', new_costs_breakdown.time_penalty);
        fprintf('  �����𻵳ɱ�: %.2f\n', new_costs_breakdown.cargo_damage);
        fprintf('  ��ʪ�ȿ��Ƴɱ�: %.2f\n', new_costs_breakdown.refrigeration);
    end
    best_costs_history(iter) = best_solution_cost; % ��¼��ǰ���������ųɱ�


    % ��Ϣ�ظ��� (����ģ�� Ant Cycle Model)
    % 1. ��Ϣ�ػӷ�
    tau = (1 - rho) * tau;

    % 2. ��Ϣ����ǿ (�������ϸ�����������������Ϣ��)
    for ant_idx = 1:num_ants
        if ~isinf(ant_costs(ant_idx)) && ant_costs(ant_idx) > 0 % ������Ч����в���
            routes_this_ant = ant_solutions{ant_idx};
            delta_tau = Q / ant_costs(ant_idx); % �ɱ�Խ�ͣ����ӵ���Ϣ��Խ��
            for r = 1:length(routes_this_ant)
                route = routes_this_ant{r};
                for i = 1:(length(route)-1)
                    node1 = route(i);
                    node2 = route(i+1);
                    tau(node1, node2) = tau(node1, node2) + delta_tau;
                    tau(node2, node1) = tau(node2, node1) + delta_tau; % �Գ�·����Ϣ����ͬ
                end
            end
        end
    end

    % (��ѡ) ��Ϣ������ (Max-Min Ant System������, ���ϸ�����ģ��)
    % tau_min = 1e-4; tau_max = 10;
    % tau(tau < tau_min) = tau_min;
    % tau(tau > tau_max) = tau_max;

end % ��������ѭ��

fprintf('ACO�㷨������\n���Ž�ɱ�: %.2f\n', best_solution_cost);
[final_cost, final_costs_breakdown, final_num_vehicles] = calculate_total_cost(best_solution_routes, customer_data, depot_data, params);
fprintf('���Ž�ĳɱ�����:\n');
fprintf('  ʹ�ó�����: %d\n', final_num_vehicles);
fprintf('  ��ʻ�ɱ�: %.2f\n', final_costs_breakdown.distance);
fprintf('  �����̶��ɱ�: %.2f\n', final_costs_breakdown.fixed_vehicle);
fprintf('  ʱ�䴰�ͷ��ɱ�: %.2f\n', final_costs_breakdown.time_penalty);
fprintf('  ����ɱ�: %.2f\n', final_costs_breakdown.cargo_damage);
fprintf('  ��ʪ�ȿ��Ƴɱ�: %.2f\n', final_costs_breakdown.refrigeration);
fprintf('  �ܳɱ� (��֤): %.2f\n', final_cost);

%% ���·�������Ŀͻ�����Ϣ
fprintf('\n����·������:\n');
if isempty(best_solution_routes)
    fprintf('δ�ҵ���Ч��·����\n');
else
    for r = 1:length(best_solution_routes)
        route = best_solution_routes{r};
        fprintf('���� %d ��·��: ', r);
        for i = 1:length(route)
            if i == length(route)
                fprintf('%d\n', route(i)); % ���һ���ڵ㲻�Ӽ�ͷ
            else
                fprintf('%d -> ', route(i)); % �ڵ�֮���ü�ͷ����
            end
        end
    end
end

%% 4. ���ӻ�
% ����ͼ: ��ʾ�㷨��������
figure;
plot(1:num_iterations, best_costs_history, 'LineWidth', 2);
xlabel('��������');
ylabel('�����ܳɱ�');
title('CVRPTW��ACO�㷨����ͼ');
grid on;

% ·��ͼ: ��ʾ���Ž�ĳ���·��
figure;
hold on;
plot(depot_data.x, depot_data.y, 'ks', 'MarkerSize', 12, 'MarkerFaceColor', 'r'); % ����Depot
text(depot_data.x+1, depot_data.y+1, 'Depot');
plot([customer_data.x], [customer_data.y], 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); % ���ƿͻ���
for i = 1:num_customers
    text(customer_data(i).x+1, customer_data(i).y+1, sprintf('%d', i)); % ��ǿͻ����
end

colors = lines(length(best_solution_routes)); % Ϊ��ͬ·�����ɲ�ͬ��ɫ
if isempty(best_solution_routes)
    disp('δ�ҵ�����������޷�����·��ͼ��');
else
    for r = 1:length(best_solution_routes) % ����ÿ������·��
        route = best_solution_routes{r};
        route_coords_x = all_nodes_coords(route, 1); % ��ȡ·���и��ڵ��x����
        route_coords_y = all_nodes_coords(route, 2); % ��ȡ·���и��ڵ��y����
        plot(route_coords_x, route_coords_y, 'LineWidth', 1.5, 'Color', colors(r,:), 'Marker', '>'); % ����·��
    end
end

axis equal; % ����x,y�����һ��
title(['�ҵ�������·�� (�ɱ�: ', sprintf('%.2f', best_solution_cost), ')']);
xlabel('X ����');
ylabel('Y ����');
legend('Depot', '�ͻ���', 'Location', 'bestoutside'); % ͼ��
hold off;

rng(s); % �ָ��������������״̬