clear; clc; close all; % ��������������ڣ��ر�����ͼ��

%% 1. ���ⶨ�� (��ACO�汾��ͬ)
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

for i=1:num_customers
    customer_data(i).late_tw = max(customer_data(i).late_tw, customer_data(i).early_tw + customer_data(i).service_time + 10);
end

depot_data.x = 35;
depot_data.y = 35;
depot_data.early_tw = 0;
depot_data.late_tw = 230; % �ֿ�����ر�ʱ��

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

%% 2. ������㷨 (NN) ʵ��
fprintf('VRPTW��������㷨��ʼ...\n');

solution_routes = {};       % �洢���ս��·��
visited_customers = false(1, num_customers); % ��ǿͻ��Ƿ��ѱ����� (�ͻ�����1��num_customers)
num_vehicles_used = 0;

while ~all(visited_customers) % �����пͻ�δ������ʱ
    num_vehicles_used = num_vehicles_used + 1;
    current_route = [1]; % ��·����Depot (�ڵ�1) ��ʼ
    current_load = 0;    % ��ǰ��������
    current_node = 1;    % ��ǰ���ڽڵ�
    % �����Ӳֿ�����ĳ�ʼʱ�䣬���ڲֿ�����ɳ���ʱ��
    departure_time_from_current_node = all_nodes_early_tw(1); 

    while true % ������������·��
        best_next_node = -1;          % �洢�ҵ��������һ���ڵ�
        min_distance_to_next = inf;   % �������һ���ڵ����С����
        
        % Ϊ�����һ���ڵ�Ԥ�����ʱ����Ϣ (�����ѡ��)
        arrival_at_best_next_cand = -1;
        service_start_at_best_next_cand = -1;
        departure_from_best_next_cand = -1;

        % Ѱ����һ���ɷ���ġ�����Ŀͻ��ڵ�
        for candidate_node_idx = 2:num_total_nodes % �������пͻ��ڵ� (����2��N+1)
            customer_actual_idx = candidate_node_idx - 1; % ��Ӧ customer_data �� visited_customers �е�����

            if ~visited_customers(customer_actual_idx) && ... % ����ÿͻ�δ��ȫ�ַ��ʹ�
               (current_load + all_nodes_demand(candidate_node_idx) <= params.vehicle_capacity) % �Ҳ�������������

                % ʱ������Լ��
                travel_time_to_candidate = dist_matrix(current_node, candidate_node_idx);
                arrival_at_candidate = departure_time_from_current_node + travel_time_to_candidate;
                
                % �絽��ȴ����������ʼʱ��
                service_start_at_candidate = max(arrival_at_candidate, all_nodes_early_tw(candidate_node_idx));
                
                % �뿪��ѡ�ͻ����ʱ��
                departure_from_candidate = service_start_at_candidate + all_nodes_service_time(candidate_node_idx);
                
                % �Ӻ�ѡ�ͻ��㷵�زֿ�����ʱ��
                time_to_return_depot_from_candidate = dist_matrix(candidate_node_idx, 1);
                
                % ��飺�������ѡ�ͻ������زֿ⣬�Ƿ�ᳬ���ֿ������ر�ʱ��
                if departure_from_candidate + time_to_return_depot_from_candidate <= all_nodes_late_tw(1) % Depot����ʱ��
                    % �������������ʱ��Լ�������Ǿ���
                    distance_to_candidate = dist_matrix(current_node, candidate_node_idx);
                    if distance_to_candidate < min_distance_to_next
                        min_distance_to_next = distance_to_candidate;
                        best_next_node = candidate_node_idx;
                        % ��¼�����ѡ�ߵ�ʱ����Ϣ����������ձ�ѡΪ���
                        arrival_at_best_next_cand = arrival_at_candidate;
                        service_start_at_best_next_cand = service_start_at_candidate;
                        departure_from_best_next_cand = departure_from_candidate;
                    end
                end
            end
        end % �����Ժ�ѡ�ڵ�ı���

        if best_next_node == -1 % ���û���ҵ����е���һ���ͻ�
            break; % ������ǰ������·������
        else
            % ���ҵ��������һ���ڵ���ӵ���ǰ·��
            current_route = [current_route, best_next_node];
            current_load = current_load + all_nodes_demand(best_next_node);
            visited_customers(best_next_node-1) = true; % ��Ǹÿͻ��ѱ����� (ע������ת��)
            current_node = best_next_node; % ���µ�ǰ�ڵ�
            
            % ���´ӵ�ǰ�ڵ㣨���¼����best_next_node��������ʱ��
            % ���ʱ���Ƿ�����best_next_node����뿪ʱ��
            departure_time_from_current_node = departure_from_best_next_cand;
        end
    end % ������������·������

    % ·��ĩβ���Depot���γ�������·
    current_route = [current_route, 1];
    if length(current_route) > 2 % ֻ�е�·������������һ���ͻ�ʱ����һ����Ч·��
        solution_routes{end+1} = current_route;
    elseif ~all(visited_customers) && num_vehicles_used > num_customers % ��ȫ�жϣ���ֹ���пͻ����޷�����ʱ��ѭ��
        fprintf('����: �޷�Ϊ����ʣ��ͻ��滮·�������ֿͻ�����δ������\n');
        break; % �˳���ѭ��
    end
    
    if isempty(solution_routes) && ~all(visited_customers) && num_customers > 0
         % ���һ��·�������������������пͻ�δ����˵������������
         % �����һ���ͻ����޷�����Լ��
         fprintf('����: δ�ܹ����κ���Ч·���������пͻ�δ����\n');
         break;
    end

end % ���� while ~all(visited_customers)

%% 3. �����ܳɱ���������
if ~all(visited_customers) && num_customers > 0
    fprintf('ע��: �������пͻ���������!\n');
    fprintf('�ѷ���ͻ���: %d/%d\n', sum(visited_customers), num_customers);
end

if isempty(solution_routes) && num_customers > 0
    fprintf('������㷨δ���ҵ��κ���Ч·����\n');
    solution_cost = inf;
    costs_breakdown = struct('distance', inf, 'fixed_vehicle', inf, 'time_penalty', inf, 'cargo_damage', inf, 'refrigeration', inf);
    actual_num_vehicles = 0;
else
    [solution_cost, costs_breakdown, actual_num_vehicles] = calculate_total_cost(solution_routes, customer_data, depot_data, params);
end

fprintf('������㷨������\n');
fprintf('�ҵ��Ľ�ĳɱ�: %.2f\n', solution_cost);
if ~isinf(solution_cost)
    fprintf('�ɱ�����:\n');
    fprintf('  ʹ�ó�����: %d\n', actual_num_vehicles);
    fprintf('  ��ʻ�ɱ�: %.2f\n', costs_breakdown.distance);
    fprintf('  �����̶��ɱ�: %.2f\n', costs_breakdown.fixed_vehicle);
    fprintf('  ʱ�䴰�ͷ��ɱ�: %.2f\n', costs_breakdown.time_penalty);
    fprintf('  ����ɱ�: %.2f\n', costs_breakdown.cargo_damage);
    fprintf('  ��ʪ�ȿ��Ƴɱ�: %.2f\n', costs_breakdown.refrigeration);
end

%% 4. ���ӻ�
% "����ͼ" ����NN���岻����Ϊ��ֻ����һ���⡣���ﻭһ�����ʾ���ճɱ���
figure;
plot(1, solution_cost, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlim([0 2]); % ����X�᷶Χ�Ա���ʾ������
xlabel('�� (������㷨)');
ylabel('�ܳɱ�');
title('������㷨���CVRPTW�ĳɱ�');
grid on;
if isinf(solution_cost)
    text(1, 0, '�޽��ɱ����޴�', 'HorizontalAlignment', 'center');
else
    text(1, solution_cost, sprintf('%.2f', solution_cost), 'VerticalAlignment','bottom', 'HorizontalAlignment','center');
end


% ·��ͼ
figure;
hold on;
plot(depot_data.x, depot_data.y, 'ks', 'MarkerSize', 12, 'MarkerFaceColor', 'r'); % ����Depot
text(depot_data.x+1, depot_data.y+1, 'Depot');
plot([customer_data.x], [customer_data.y], 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); % ���ƿͻ���
for i = 1:num_customers
    text(customer_data(i).x+1, customer_data(i).y+1, sprintf('%d', i)); % ��ǿͻ����
end

if isempty(solution_routes) || isinf(solution_cost)
    disp('����Ч·���ɹ����ơ�');
else
    colors = lines(length(solution_routes)); % Ϊ��ͬ·�����ɲ�ͬ��ɫ
    for r = 1:length(solution_routes) % ����ÿ��·��
        route = solution_routes{r};
        route_coords_x = all_nodes_coords(route, 1); % ��ȡ·���и��ڵ��x����
        route_coords_y = all_nodes_coords(route, 2); % ��ȡ·���и��ڵ��y����
        plot(route_coords_x, route_coords_y, 'LineWidth', 1.5, 'Color', colors(r,:), 'Marker', '>'); % ����·��
    end
end

axis equal; % ����x,y�����һ��
if isinf(solution_cost)
    title('������㷨·��ͼ (����Ч��)');
else
    title(['������㷨�ҵ���·�� (�ɱ�: ', sprintf('%.2f', solution_cost), ')']);
end
xlabel('X ����');
ylabel('Y ����');
legend('Depot', '�ͻ���', 'Location', 'bestoutside'); % ͼ��
hold off;

% ��·���滮��ɺ����ÿ������·��
if isempty(solution_routes) || isinf(solution_cost)
    fprintf('δ�ҵ���Ч·����\n');
else
    fprintf('���г�����·���滮���£�\n');
    for r = 1:length(solution_routes)
        route = solution_routes{r};
        fprintf('���� %d ��·��: ', r);
        % ��������ӡ·���ϵĽڵ���
        for i = 1:length(route)
            if i == length(route)
                fprintf('%d\n', route(i)); % ���һ���ڵ㻻��
            else
                fprintf('%d -> ', route(i)); % ���ӷ���
            end
        end
    end
end

rng(s); % �ָ��������������״̬