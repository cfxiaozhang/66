% main_cvrptw_aco.m

% --- 0. ��������������������� ---
clear; clc; close all;
rng('default'); % Ϊ�˽���ɸ���

% --- 1. ���ⶨ�� ---
% �ڵ�1�ǲֿ� (Depot)���ڵ�2��N+1�ǿͻ���
% Coords = [x_depot, y_depot; x_cust1, y_cust1; ...];
% Demand = [0; demand_cust1; ...]; % �ֿ�����Ϊ0
% ServiceTime = [0; service_cust1; ...]; % �ֿ����ʱ��Ϊ0
% TimeWindow = [ES_depot, LS_depot; ES_cust1, LS_cust1; ...]; % [���翪ʼ����ʱ��, ����ʼ����ʱ��]

% ʾ������ (���滻Ϊ����ʵ������)
num_customers = 20; % �ͻ�����
depot_coord = [50 50]; % �ֿ�����
customer_coords = rand(num_customers, 2) * 100; % �ͻ����� (0-100��Χ)

Nodes.Coords = [depot_coord; customer_coords];
Nodes.Demand = [0; randi([10, 30], num_customers, 1)]; % �ͻ�����
Nodes.ServiceTime = [0; randi([5, 10], num_customers, 1)]; % �ͻ�����ʱ��

% ʱ�䴰����: ES_i = ���, LS_i = ES_i + ���ԣ��
Nodes.TimeWindow = [0, 1000; % �ֿ��ʱ�䴰 (ͨ���ܴ�)
                   randi([0, 200], num_customers, 1), zeros(num_customers, 1)]; 
for i = 1:num_customers
    Nodes.TimeWindow(i+1, 2) = Nodes.TimeWindow(i+1, 1) + randi([50, 150]); % LS = ES + ԣ��
end

Vehicle.Capacity = 150; % ��������
Vehicle.Speed = 1;     % �����ٶ� (���뵥λ/ʱ�䵥λ)

% �ɱ�����
Cost.PerDistance = 1.0;     % ÿ��λ����ĳɱ�
Cost.PerWaitTime = 0.5;     % ÿ��λ�ȴ�ʱ��ĳɱ�
Cost.PerTardinessTime = 2.0;% ÿ��λ��ʱ��ĳͷ��ɱ� (����ʼʱ������LS)

num_nodes = size(Nodes.Coords, 1); % �ܽڵ��� (�ֿ� + �ͻ�)

% --- 2. ACO ���� ---
ACO.NumAnts = 30;        % ��������
ACO.MaxIter = 150;       % ����������
ACO.Alpha = 1;           % ��Ϣ���������� (Pheromone influence)
ACO.Beta = 3;            % ������������ (Heuristic influence)
ACO.Rho = 0.1;           % ��Ϣ�������� (Pheromone evaporation rate)
ACO.Q = 100;             % ��Ϣ��ǿ�ȳ��� (Pheromone deposit constant)
ACO.Tau0 = 0.1;          % ��ʼ��Ϣ��ˮƽ (Initial pheromone level)

% --- 3. ��ʼ�� ---
% ����������
DistMatrix = zeros(num_nodes, num_nodes);
for i = 1:num_nodes
    for j = i+1:num_nodes
        DistMatrix(i,j) = norm(Nodes.Coords(i,:) - Nodes.Coords(j,:));
        DistMatrix(j,i) = DistMatrix(i,j);
    end
end

% ��ʼ����Ϣ�ؾ���
Pheromone = ones(num_nodes, num_nodes) * ACO.Tau0;

% ����������Ϣ���� (ͨ���Ǿ���ĵ���)
Heuristic = zeros(num_nodes, num_nodes);
for i = 1:num_nodes
    for j = 1:num_nodes
        if DistMatrix(i,j) > 0
            Heuristic(i,j) = 1 / DistMatrix(i,j); 
        end
    end
end

% �洢���Ž�
BestSol.Routes = {};
BestSol.Cost = inf;
BestSol.CostHistory = zeros(ACO.MaxIter, 1); % ��¼ÿ�ε��������ųɱ�

% --- 4. ��ACOѭ�� ---
fprintf('��ʼACO���CVRPTW...\n');
for iter = 1:ACO.MaxIter
    AntSolutions = cell(ACO.NumAnts, 1); % �洢ÿֻ���ϵĽ�
    AntCosts = zeros(ACO.NumAnts, 1);    % �洢ÿֻ���Ͻ���ܳɱ�

    for k = 1:ACO.NumAnts % ��ÿֻ����
        CurrentAnt.Routes = {};    % ��ǰ���Ϲ�����·������
        CurrentAnt.TotalCost = 0;  % ��ǰ���ϵ��ܳɱ�
        
        CustomersToVisit = true(1, num_customers); % ��ǿͻ��Ƿ��ѱ����� (true��ʾδ����)
        num_unvisited = num_customers; % δ���ʿͻ�����

        while num_unvisited > 0 % �����пͻ�δ������ʱ
            % Ϊ��ǰ���Ͽ�ʼһ����·�� (һ���³�)
            current_route = [1]; % �Ӳֿ���� (�ڵ�����1)
            current_load = 0;    % ��ǰ��������
            current_time = Nodes.TimeWindow(1,1); % �����Ӳֿ������ʱ�� (ͨ���ǲֿ��ES, e.g., 0)
            current_route_cost = 0; % ��ǰ·���ĳɱ�

            % ������ǰ·��
            while true 
                last_node_idx = current_route(end); % ��ǰ·���е����һ���ڵ�
                
                possible_next_nodes_indices = []; % �洢��ѡ��һ���ڵ��ʵ������
                selection_probabilities = [];     % �洢ѡ�����ѡ�ڵ�ĸ���

                % �������пͻ���Ѱ�ҿ��е���һ���ͻ�
                for cust_idx = 1:num_customers % cust_idx �ǿͻ��ı�� (1 �� num_customers)
                    actual_node_idx = cust_idx + 1; % actual_node_idx �ǿͻ���Nodes�ṹ�е�����

                    if CustomersToVisit(cust_idx) % ����ͻ���δ������
                        % 1. �������Լ��
                        if current_load + Nodes.Demand(actual_node_idx) <= Vehicle.Capacity
                            % 2. ����ѡ��ýڵ��������Ϣ����Ϣ��Ũ��
                            % (��ʱ�䴰��ζ������"����"�����ɱ���仯���˴�ѡ����ʲ�ֱ�ӳͷ����ɱ���֮�����)
                            prob_value = (Pheromone(last_node_idx, actual_node_idx) ^ ACO.Alpha) * ...
                                         (Heuristic(last_node_idx, actual_node_idx) ^ ACO.Beta);
                            
                            if prob_value == 0 || isnan(prob_value) || isinf(prob_value)
                                prob_value = 1e-9; % ��ֹ����Ϊ0����Чֵ
                            end
                            
                            possible_next_nodes_indices(end+1) = actual_node_idx;
                            selection_probabilities(end+1) = prob_value;
                        end
                    end
                end

                if isempty(possible_next_nodes_indices)
                    break; % û�п�ѡ�����һ���ͻ� (�����пͻ��ѷ������)��������ǰ·������
                end

                % ���̶�ѡ����һ���ͻ�
                selection_probabilities = selection_probabilities / sum(selection_probabilities);
                if any(isnan(selection_probabilities)) || any(isinf(selection_probabilities)) || sum(selection_probabilities)==0
                    rand_idx = randi(length(possible_next_nodes_indices)); % ���ѡһ��
                else
                    r = rand();
                    cum_probs = cumsum(selection_probabilities);
                    chosen_idx_in_possible = find(r <= cum_probs, 1, 'first');
                    if isempty(chosen_idx_in_possible) % �Է���һ
                         rand_idx = randi(length(possible_next_nodes_indices));
                         chosen_idx_in_possible = rand_idx;
                    end
                end
                chosen_node_actual_idx = possible_next_nodes_indices(chosen_idx_in_possible);
                
                % --- ��ѡ�еĿͻ����뵱ǰ·����������·��״̬ ---
                % a. ������ʻ�ɱ���ʱ��
                travel_dist = DistMatrix(last_node_idx, chosen_node_actual_idx);
                travel_time = travel_dist / Vehicle.Speed;
                current_route_cost = current_route_cost + travel_dist * Cost.PerDistance;
                
                arrival_at_chosen = current_time + travel_time; % ����ѡ�пͻ���ʱ��
                
                % b. �������ʱ�䡢�ȴ��ɱ������ɱ�
                service_start_chosen = max(arrival_at_chosen, Nodes.TimeWindow(chosen_node_actual_idx, 1));
                
                wait_time_chosen = service_start_chosen - arrival_at_chosen;
                current_route_cost = current_route_cost + wait_time_chosen * Cost.PerWaitTime;
                
                tardiness_chosen = max(0, service_start_chosen - Nodes.TimeWindow(chosen_node_actual_idx, 2));
                current_route_cost = current_route_cost + tardiness_chosen * Cost.PerTardinessTime;
                
                departure_time_chosen = service_start_chosen + Nodes.ServiceTime(chosen_node_actual_idx);
                
                % c. ����·��״̬
                current_route(end+1) = chosen_node_actual_idx; % ��ӵ�·��
                current_load = current_load + Nodes.Demand(chosen_node_actual_idx); % ��������
                current_time = departure_time_chosen; % ���µ�ǰʱ��Ϊ�뿪�˿ͻ���ʱ��
                
                % d. ��ǿͻ�Ϊ�ѷ���
                CustomersToVisit(chosen_node_actual_idx - 1) = false; % -1 ����ΪCustomersToVisit������Ӧ�ͻ����
                num_unvisited = num_unvisited - 1;
            end % ������ǰ·���Ŀͻ���� (inner while true)
            
            % ��ǰ·��������� (û�и���ͻ��ɷ��� �� �������������Է����κ�ʣ��ͻ�)
            % ���·������������һ���ͻ����򷵻زֿ�
            if length(current_route) > 1 
                last_customer_in_route_idx = current_route(end);
                
                % a. ���㷵�زֿ����ʻ�ɱ���ʱ��
                travel_dist_to_depot = DistMatrix(last_customer_in_route_idx, 1);
                travel_time_to_depot = travel_dist_to_depot / Vehicle.Speed;
                current_route_cost = current_route_cost + travel_dist_to_depot * Cost.PerDistance;
                
                arrival_at_depot = current_time + travel_time_to_depot;
                
                % b. ��鷵�زֿ��Ƿ���� (����ֿ���������ʱ��LS)
                depot_tardiness = max(0, arrival_at_depot - Nodes.TimeWindow(1,2));
                current_route_cost = current_route_cost + depot_tardiness * Cost.PerTardinessTime;
                
                current_route(end+1) = 1; % ·��ĩβ��Ӳֿ�
                
                CurrentAnt.Routes{end+1} = current_route; % ����·���������ϵĽ���
                CurrentAnt.TotalCost = CurrentAnt.TotalCost + current_route_cost; %�ۼ��ܳɱ�
            elseif num_unvisited > 0 && isempty(CurrentAnt.Routes) && length(current_route) == 1
                % �������: �����޷������κοͻ� (�����ǲ��������ⶨ������)
                % ������ʱ�䴰����Ӧ�ú��ٷ�������������Զ��������ʱ�䴰���Ȳ�����
                CurrentAnt.TotalCost = inf; % ���輫�߳ɱ�
                break; % �˳���ǰ���ϵ�·������
            end

        end % ������ǰ���ϵ�����·������ (while num_unvisited > 0)
        
        AntSolutions{k} = CurrentAnt; % �洢��ǰ���ϵ�������
        AntCosts(k) = CurrentAnt.TotalCost; % �洢��ǰ���Ͻ���ܳɱ�

        % ����ȫ�����Ž�
        if CurrentAnt.TotalCost < BestSol.Cost
            BestSol.Cost = CurrentAnt.TotalCost;
            BestSol.Routes = CurrentAnt.Routes;
            fprintf('���� %d, ���� %d: �����ųɱ� = %.2f\n', iter, k, BestSol.Cost);
        end
    end % ��������ѭ��

    % --- ��Ϣ�ظ��� ---
    % 1. ��Ϣ������
    Pheromone = (1 - ACO.Rho) * Pheromone;

    % 2. ��Ϣ�س���
    for k_ant = 1:ACO.NumAnts
        if isinf(AntCosts(k_ant)); continue; end % ������Ͻ���Ч��������
        
        ant_solution_routes = AntSolutions{k_ant}.Routes;
        delta_tau = ACO.Q / AntCosts(k_ant); % Q / Lk
        
        for r = 1:length(ant_solution_routes) % �������ϵ�ÿ��·��
            route = ant_solution_routes{r};
            for i = 1:length(route)-1 % ����·���ϵ�ÿ����
                node_from = route(i);
                node_to = route(i+1);
                Pheromone(node_from, node_to) = Pheromone(node_from, node_to) + delta_tau;
                Pheromone(node_to, node_from) = Pheromone(node_to, node_from) + delta_tau; % �Գ�����
            end
        end
    end
    
    % ȷ����Ϣ�ز�����ĳ����ֵ����ֹ����ͣ��
    Pheromone(Pheromone < ACO.Tau0/10) = ACO.Tau0/10; % ���磬��ʼ��Ϣ�ص�1/10

    BestSol.CostHistory(iter) = BestSol.Cost; % ��¼���ε��������ųɱ�
    if mod(iter, 10) == 0 || iter == 1 || iter == ACO.MaxIter
      fprintf('���� %d: ��ǰ���ųɱ� = %.2f\n', iter, BestSol.Cost);
    end
end % ��������ѭ��

fprintf('ACO�����ɡ�\n�ҵ������ųɱ�: %.2f\n', BestSol.Cost);

% --- 5. �������ӻ� ---

% ����ͼ
figure;
plot(1:ACO.MaxIter, BestSol.CostHistory, 'LineWidth', 2);
xlabel('��������');
ylabel('�ܳɱ�');
title('ACO�������� (CVRPTW)');
grid on;

% ·��ͼ
figure;
hold on;
% ���Ʋֿ�
plot(Nodes.Coords(1,1), Nodes.Coords(1,2), 'ks', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'DisplayName', '�ֿ�');
% ���ƿͻ���
plot(Nodes.Coords(2:end,1), Nodes.Coords(2:end,2), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', '�ͻ�');

% Ϊÿ��·��ʹ�ò�ͬ��ɫ
route_colors = lines(length(BestSol.Routes)); 
for r = 1:length(BestSol.Routes)
    route_nodes_indices = BestSol.Routes{r}; % ��ȡ·���еĽڵ�����
    route_coords = Nodes.Coords(route_nodes_indices, :); % ��ȡ��Щ�ڵ������
    plot(route_coords(:,1), route_coords(:,2), 'LineWidth', 1.5, 'Color', route_colors(r,:), 'Marker', 'none');
    
    % Ϊ�ͻ��ڵ�����ı���ǩ (�ͻ����, ����, ʱ�䴰)
    for i = 1:length(route_nodes_indices)
        node_idx = route_nodes_indices(i);
        if node_idx > 1 % ����ǿͻ��ڵ�
            customer_label_num = node_idx - 1;
            text_str = sprintf('C%d (D:%d, TW:[%d,%d])', ...
                               customer_label_num, ...
                               Nodes.Demand(node_idx), ...
                               Nodes.TimeWindow(node_idx,1), ...
                               Nodes.TimeWindow(node_idx,2));
            text(Nodes.Coords(node_idx,1)+0.8, Nodes.Coords(node_idx,2)+0.8, text_str, 'FontSize', 7);
        end
    end
end

title(['�ҵ�������·������ (�ܳɱ�: ' num2str(BestSol.Cost, '%.2f') ')']);
xlabel('X ����');
ylabel('Y ����');
legend('show', 'Location', 'bestoutside');
axis equal;
grid on;
hold off;

% �ڿ���̨�������·��
disp('�ҵ�������·��:');
for r = 1:length(BestSol.Routes)
    % ���ڵ�����ӳ��ؿͻ���� (0 ����ֿ�)
    route_display = BestSol.Routes{r} - 1; 
    fprintf('·�� %d: %s\n', r, mat2str(route_display));
end

% (��ѡ) ��ϸ��֤���Ž�ĳɱ�����
TotalVerifiedCost = 0;
disp('���Ž����ϸ�ɱ��ֽ�:');
for r_idx = 1:length(BestSol.Routes)
    route = BestSol.Routes{r_idx};
    route_cost_detail = 0;
    current_time_verify = Nodes.TimeWindow(1,1); % �Ӳֿ����ʱ��
    current_load_verify = 0;
    
    fprintf('\n·�� %d: �ֿ� -> ', r_idx);
    
    for i = 1:(length(route)-1) % ����·����ÿһ��
        from_node = route(i);
        to_node = route(i+1);
        
        % 1. ��ʻ�ɱ�
        dist = DistMatrix(from_node, to_node);
        travel_c = dist * Cost.PerDistance;
        route_cost_detail = route_cost_detail + travel_c;
        
        travel_t = dist / Vehicle.Speed;
        arrival_at_to = current_time_verify + travel_t;
        
        if to_node ~= 1 % ���Ŀ���ǿͻ�
            fprintf('C%d -> ', to_node-1);
            
            % 2. �ȴ��ɱ�
            service_start = max(arrival_at_to, Nodes.TimeWindow(to_node,1));
            wait_t = service_start - arrival_at_to;
            wait_c = wait_t * Cost.PerWaitTime;
            route_cost_detail = route_cost_detail + wait_c;
            
            % 3. ���ͷ��ɱ�
            tardiness_t = max(0, service_start - Nodes.TimeWindow(to_node,2));
            tardiness_c = tardiness_t * Cost.PerTardinessTime;
            route_cost_detail = route_cost_detail + tardiness_c;
            
            % ����ʱ��͸���
            departure_time = service_start + Nodes.ServiceTime(to_node);
            current_time_verify = departure_time;
            current_load_verify = current_load_verify + Nodes.Demand(to_node);

            fprintf('(����:%.1f, �ȴ�:%.1f, ����ʼ:%.1f, ��:%.1f, �뿪:%.1f, ����:%d) - ', ...
                arrival_at_to, wait_t, service_start, tardiness_t, departure_time, current_load_verify);

        else % ���Ŀ���ǲֿ� (·������)
            fprintf('�ֿ� ');
            current_time_verify = arrival_at_to; % ����ֿ�ʱ��
            % ��鷵�زֿ��Ƿ����
            depot_tardiness_t = max(0, current_time_verify - Nodes.TimeWindow(1,2));
            depot_tardiness_c = depot_tardiness_t * Cost.PerTardinessTime;
            route_cost_detail = route_cost_detail + depot_tardiness_c;
            fprintf('(����ֿ�:%.1f, �ֿ���:%.1f) ', current_time_verify, depot_tardiness_t);
        end
    end
    fprintf('\n·�� %d �ɱ�: %.2f. ��������: %d (��������: %d)\n', r_idx, route_cost_detail, current_load_verify, Vehicle.Capacity);
    TotalVerifiedCost = TotalVerifiedCost + route_cost_detail;
end
fprintf('\n����֤�ɱ� (������ϸ�ֽ�): %.2f\n', TotalVerifiedCost);
fprintf('ACO�ҵ������ųɱ�: %.2f\n', BestSol.Cost);
if abs(TotalVerifiedCost - BestSol.Cost) < 1e-3
    fprintf('�ɱ���֤�ɹ���\n');
else
    fprintf('����: �ɱ���֤ʧ�ܡ����: %.4f\n', abs(TotalVerifiedCost - BestSol.Cost));
end