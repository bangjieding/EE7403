var_table = load('./experiments.mat')


figure
plot(t, x,'red--');  % 绘制正弦曲线
hold on;  % 将正弦曲线保持在图形中
plot(t, y,'y+'); % 绘制余弦曲线，完成后图形中就会同时显示正弦曲线和余弦曲线
hold on
z=x+y;
plot(t,z,'b:.')
hold on
g=x.*y;
plot(t,g,'g.')
hold on
legend('sin(t)','cos(t)','sin(t)+cos(t)','sin(t)*cos(t)') %可依次设置成你想要的名字


% tem = 3:2:13;
% for i = 1:6
% 	[min_mean_val(i), min_mean_idx(i)] = min(tuning(:, tem(i)));
% 	% [min_stdvar_val(i), min_stdvar_idx(i)] = min(tuning(:, tem(i)+1));
% end
% ≈
% min_mean_line = tuning(min_mean_idx, :);
% min_stdvar_line = tuning(min_stdvar_idx, :);

% min_stdvar_line_2 = tuning(min_stdvar_idx, :);
% for i = 1:6
% 	[min_mean_val(i), min_mean_idx(i)] = min(tuning(:, tem(i)));
% 	% [min_stdvar_val(i), min_stdvar_idx(i)] = min(tuning(:, tem(i)+1));
% end


% min_mean_idx = min_mean_idx';
% min_stdvar_idx = min_stdvar_idx';