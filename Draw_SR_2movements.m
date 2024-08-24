% Authors: Lin, Li
% control variables:
p1 = [0.3,0.7]; % probability % p1 and s1 need to have the same length; % element of p1 needs to be >0
p2 = [0.5,0.5]; % probability % p2 and s2 need to have the same length; % element of p2 needs to be >0
s1 = [1,2]; % I-SFR % element of s1 needs to be >=0
s2 = [1,2]; % I-SFR % element of s2 needs to be >=0
theta = [0.5;0.5];%prediction ability of [movement1, movement2] % 0 <= theta_i <= 1
% [1;0] indicates O={1}, [0;1] indicates O={2}; [1;1] indicates O={1,2},and [0;0] indicates O = \emptyset
K_max = 1; % 1-K_max equals the ratio of all-red
%---------------------------
Total_split = 101; 
Pe = nan(1,length(p1)*length(p2)); 
nP = length(Pe);
A = [eye(nP),eye(nP)];
b = repmat(K_max,nP,1);
s1e = Pe; 
s2e = Pe; 
count_p = 0;
for i = 1:length(p1)
    for j = 1:length(p2)
        count_p = count_p + 1;
        Pe(count_p) = p1(i)*p2(j);
        s1e(count_p) = s1(i);
        s2e(count_p) = s2(j);
    end
end

% 注，如果不独立，Pe, s1e, s2e 需要独立输入
Pe = Pe/sum(Pe);
%s1e = [1,1,2,2];
%s2e = [1,2,1,2];
s1_ba = Pe*s1e';
s2_ba = Pe*s2e';

pre_target = [Pe.*s1e,Pe.*s2e]; 
piece_1 = 1:(length(pre_target)/2); %两部分
piece_2 = (length(pre_target)/2+1):length(pre_target);
theta = [theta(1+mod(0,length(theta))),theta(1+mod(1,length(theta)))];

pre_target_non = [Pe*s1_ba,Pe.*s2_ba];
% pre_target_theta = pre_target*theta + pre_target_non*(1-theta);
pre_target_theta = [Pe.*s1e*theta(1),Pe.*s2e*theta(2)]+...
    + [Pe.*s1_ba*(1-theta(1)),Pe.*s2_ba*(1-theta(2))];
%------------------------------------- 1. 原始问题 （全知问题）
Lamda1 = nan(Total_split,2); % 最大流量
fval_mat1 = nan(Total_split,1); 

for i=1:Total_split
    r = (i-1)/(Total_split-1);
    target = [pre_target(piece_1)*r,pre_target(piece_2)*(1-r)]; % coefficients when max the object 
    [X,fval] = linprog(-target,A,b,[],[],zeros(1,2*nP),[]);
    Lamda1(i,1) = pre_target(piece_1)*X(piece_1);
    Lamda1(i,2) = pre_target(piece_2)*X(piece_2);
    fval_mat1(i) = fval;
end
%------------------------------------- 2. 无知问题
Lamda2 = nan(Total_split,2); % 最大流量
fval_mat2 = nan(Total_split,1); 

for i=1:Total_split
    r = (i-1)/(Total_split-1);
    target = [pre_target_non(piece_1)*r,pre_target_non(piece_2)*(1-r)]; % coefficients when max the object 
    [X,fval] = linprog(-target,A,b,[],[],zeros(1,2*nP),[]);
    Lamda2(i,1) = pre_target_non(piece_1)*X(piece_1);
    Lamda2(i,2) = pre_target_non(piece_2)*X(piece_2);
    fval_mat2(i) = fval;
end
%------------------------------------- 3. 半知问题 
Lamda3 = nan(Total_split,2); % 最大流量
fval_mat3 = nan(Total_split,1); 

for i=1:Total_split
    r = (i-1)/(Total_split-1);
    target = [pre_target_theta(piece_1)*r,pre_target_theta(piece_2)*(1-r)]; % coefficients when max the object 
    [X,fval] = linprog(-target,A,b,[],[],zeros(1,2*nP),[]);
    Lamda3(i,1) = pre_target_theta(piece_1)*X(piece_1);
    Lamda3(i,2) = pre_target_theta(piece_2)*X(piece_2);
    fval_mat3(i) = fval;
end


c1 = unique(Lamda1,'rows');
c2 = unique(Lamda2,'rows');
c3 = unique(Lamda3,'rows');
scatters1 = sortrows(c1,[-1, 2]);
scatters2 = sortrows(c2,[-1, 2]);
scatters3 = sortrows(c3,[-1, 2]);
plot(scatters1(:,1),scatters1(:,2),LineWidth= 1.5);
hold on;
plot(scatters3(:,1),scatters3(:,2),LineWidth= 1.5);
hold on;
plot(scatters2(:,1),scatters2(:,2),LineWidth= 1.5);
hold off;

legend('Complete info.','Partial info. with \theta = ['+string(theta(1))'+';'+string(theta(2))+']',' No info.')
xlabel('Capacity of movement 1: veh/interval','rot', 0)
ylabel('Capacity of movement 2: veh/interval','rot', 90)
set(gca,'FontSize',14,'Fontname','Times New Roman')
% xLimits = xlim; % 获取当前x轴的范围  
% xticks(xLimits(1):0.5:xLimits(2)); % 从最小值到最大值，间隔为0.5  
