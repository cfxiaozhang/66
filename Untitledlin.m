%% ������㷨������ʱ�䴰�ĳ���·������(VRPTW)
clear;
clc;
close all;

%% �����������
n = 50;                 % �ͻ�����
m = 25;                  % ��������
capacity = 200;          % ��������
depot = [35, 35];         % ������������
F = 300;                 % �����̶�ʹ�óɱ�
C = 2;                  % ��λ������ʻ�ɱ�
mu = 0.001;               % ��λʱ�����ϵ��
P = 10000;                 % ��λ���̼�ֵ
G = 0.5;                  % ��ʪ���豸��λʱ��ɱ�
late_penalty = 2;       % ���ͷ�ϵ��
speed = 60;              % ������ʻ�ٶ�

%% ��������ͻ�����
rng(0); % ����������ӣ���֤���������
customers.x = xlsread('�ͻ���.xlsx','B11:B60');      % �ͻ�x����
customers.y = xlsread('�ͻ���.xlsx','C11:C60');      % �ͻ�y����
customers.demand = xlsread('�ͻ���.xlsx','D11:D60'); % �ͻ�����
customers.ready_time = xlsread('�ͻ���.xlsx','E11:E60'); % �ͻ��������ʱ��
customers.due_time = xlsread('�ͻ���.xlsx','F11:F60'); % �ͻ��������ʱ��
customers.service_time = xlsread('�ͻ���.xlsx','G11:G60'); % �ͻ�����ʱ��

%% ����������
distance = zeros(n+1, n+1); % ������������
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

%% ������㷨���
unvisited = 2:n+1; % δ���ʵĿͻ���2��n+1��1���������ģ�
routes = cell(m, 1); % m������·��
vehicle_load = zeros(m, 1); % ÿ�����ĵ�ǰ����
vehicle_time = zeros(m, 1); % ÿ�����ĵ�ǰʱ��
total_cost = 0; % �ܳɱ�

for k = 1:m
    current_node = 1; % ���������Ŀ�ʼ
    routes{k} = [current_node];
    
    while ~isempty(unvisited)
        nearest_node = 0;
        min_distance = inf;
        
        % Ѱ������Ŀ��пͻ�
        for i = 1:length(unvisited)
            j = unvisited(i);
            
            % �������Լ��
            if vehicle_load(k) + customers.demand(j-1) > capacity
                continue;
            end
            
            % ���㵽��ʱ���ʱ�䴰Լ��
            travel_time = distance(current_node, j) / speed;
            arrival_time = vehicle_time(k) + travel_time;
            
            % �������ͷ�����ʱ�䴰��
            late_penalty_cost = late_penalty * max(0, arrival_time - customers.due_time(j-1));
            
            % �����ܳɱ�������+���ͷ���
            total_cost_candidate = distance(current_node, j) + late_penalty_cost;
            
            % �������������Լ������������ڵ�
            if total_cost_candidate < min_distance
                min_distance = total_cost_candidate;
                nearest_node = j;
            end
        end
        
        % ���û���ҵ����пͻ���������ǰ������·��
        if nearest_node == 0
            break;
        end
        
        % ���³���״̬
        vehicle_load(k) = vehicle_load(k) + customers.demand(nearest_node-1);
        travel_time = distance(current_node, nearest_node) / speed;
        vehicle_time(k) = vehicle_time(k) + travel_time;
        wait_time = max(0, customers.ready_time(nearest_node-1) - vehicle_time(k));
        vehicle_time(k) = vehicle_time(k) + wait_time + customers.service_time(nearest_node-1);
        
        % ����·��
        routes{k} = [routes{k}, nearest_node];
        
        % ��δ���ʿͻ����Ƴ�
        unvisited(find(unvisited == nearest_node)) = [];
        
        % ���µ�ǰ�ڵ�
        current_node = nearest_node;
    end
    
    % ����������������
    if length(routes{k}) > 1 % ��������з���ͻ�
        routes{k} = [routes{k}, 1];
        
        % ���㳵�������������ĵ�ʱ��
        travel_time = distance(current_node, 1) / speed;
        vehicle_time(k) = vehicle_time(k) + travel_time;
    end
end

%% �����ܳɱ�
total_cost = 0;
total_distance = 0;
total_fixed_cost = 0;
total_late_penalty = 0;
total_damage_cost = 0;
total_equipment_cost = 0;

for k = 1:m
    route = routes{k};
    if length(route) > 2 % �������������ĺ�һ���ͻ�
        vehicle_time = 0;
        vehicle_distance = 0;
        damage_cost = 0;
        equipment_cost = 0;
        late_cost = 0;
        
        total_fixed_cost = total_fixed_cost + F;
        
        for i = 1:length(route)-1
            current_node = route(i);
            next_node = route(i+1);
            
            % ��ʻ�ɱ�
            vehicle_distance = vehicle_distance + distance(current_node, next_node);
            
            % ����ʱ�䣨���ڻ�����豸�ɱ���
            travel_time = distance(current_node, next_node) / speed;
            equipment_cost = equipment_cost + travel_time * G;
            
            if next_node ~= 1 % ���������Ľڵ㣨�ͻ��ڵ㣩
                j = next_node - 1; % �ͻ�����
                damage_cost = damage_cost + travel_time * mu * P * customers.demand(j);
                
                % ʱ�䴰�ͷ��������ͷ���
                vehicle_time = vehicle_time + travel_time;
                arrival_time = vehicle_time;
                late_cost = late_cost + late_penalty * max(0, arrival_time - customers.due_time(j));
                
                % ���³���ʱ��
                wait_time = max(0, customers.ready_time(j) - arrival_time);
                vehicle_time = arrival_time + wait_time + customers.service_time(j);
            else
                vehicle_time = vehicle_time + travel_time; % �����������ĵ�ʱ��
            end
        end
        
        total_distance = total_distance + vehicle_distance;
        total_late_penalty = total_late_penalty + late_cost;
        total_damage_cost = total_damage_cost + damage_cost;
        total_equipment_cost = total_equipment_cost + equipment_cost;
    end
end

% �����ܳɱ�
total_cost = total_distance * C + total_fixed_cost + total_late_penalty + total_damage_cost + total_equipment_cost;

%% ������
fprintf('������㷨�����:\n');
fprintf('�ܳɱ�: %.2f\n', total_cost);

% ��ʾ����·��
fprintf('����·��:\n');
for k = 1:m
    route = routes{k};
    if length(route) > 2 % �������������ĺ�һ���ͻ�
        fprintf('���� %d: ', k);
        for i = 1:length(route)
            if route(i) == 1
                fprintf('��������');
            else
                fprintf('�ͻ� %d', route(i)-1);
            end
            if i < length(route)
                fprintf(' -> ');
            end
        end
        fprintf('\n');
    end
end

% �ֽ���ʾ����ɱ�
fprintf('\n�ɱ��ֽ�:\n');
fprintf('��ʻ�ɱ� (�����C): %.2f\n', total_distance * C);
fprintf('�����̶��ɱ�: %.2f\n', total_fixed_cost);
fprintf('���ͷ��ɱ�: %.2f\n', total_late_penalty);
fprintf('����ɱ�: %.2f\n', total_damage_cost);
fprintf('��ʪ���豸�ɱ�: %.2f\n', total_equipment_cost);
fprintf('�ܳɱ�: %.2f\n', total_cost);

%% ���ӻ�
% ����·��ͼ
figure;
hold on;
% ������������
plot(depot(1), depot(2), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
text(depot(1)+1, depot(2)+1, '��������');

% ���ƿͻ���
for i = 1:n
    plot(customers.x(i), customers.y(i), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    text(customers.x(i)+1, customers.y(i)+1, num2str(i));
end

% ����·��
colors = lines(m);
for k = 1:m
    route = routes{k};
    if length(route) > 2 % �������������ĺ�һ���ͻ�
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
        text(mean(x), mean(y), ['���� ', num2str(k)], 'Color', colors(k,:));
    end
end

title('������㷨���ų���·������');
xlabel('X����');
ylabel('Y����');
grid on;
axis equal;
legend('��������', '�ͻ�', '����·��');
hold off;    