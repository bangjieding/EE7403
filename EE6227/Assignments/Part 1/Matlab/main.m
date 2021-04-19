% clear all
% mex cec17_func.cpp -DWINDOWS

fhd = str2func('cec17_func');

func_num = 1;
D = 10;
Xmin = -100;
Xmax = 100;
pop_size = 100;
iter_max = 5000;
Max_FES = pop_size * iter_max - 1 ;
% C_set = [1.7, 1.3];

pso_c_set = [1.7, 1.3];
pso_w_set = [0.9, 0.3];
pso_Inertia_Weight = pso_w_set(1) - (1:iter_max) .* (pso_w_set(2) ./ iter_max);

clpso_c_set = [1.5, 1.5];
clpso_w_set = [0.8, 0.3];
clpso_Inertia_Weight = clpso_w_set(1) - (1:iter_max) .* (clpso_w_set(2) ./ iter_max);


experiment_fun_sele = [3, 4, 5, 6, 7, 8, 9, 10, 11, 21, 22];

experiment_run = 5;

for i = 1:length(experiment_fun_sele)
	func_num = experiment_fun_sele(i);
	for j = 1:experiment_run
		i, j
        [pso_gbest, pso_gbestval, pso_FES, pso_convergence] =  PSO_func(fhd, D, pop_size, iter_max, Xmin, Xmax, pso_Inertia_Weight, pso_c_set, func_num);

		[clpso_gbest, clpso_gbestval, clpso_FES, clpso_convergence] =  CLPSO_func(fhd, D, pop_size, iter_max, Max_FES, Xmin, Xmax, clpso_Inertia_Weight, clpso_c_set, func_num);

		pso_xbest(j, :) = pso_gbest;
		pso_fbest(i, j) = pso_gbestval;

		clpso_xbest(j, :) = clpso_gbest;
		clpso_fbest(i, j) = clpso_gbestval;
	end
	pso_f_mean(i) = mean(pso_fbest(i, :));
	pso_f_stdvar(i) = std(pso_fbest(i, :));

	clpso_f_mean(i) = mean(clpso_fbest(i, :));
	clpso_f_stdvar(i) = std(clpso_fbest(i, :));

	pso_con(i, :) = pso_convergence;
	clpso_con(i, :) = clpso_convergence;
end

save('experiments.mat');


% runs = 1;


% C1 = [1:0.1:2];
% % C2 = [1:0.1:2];
% % wmax = [0.1:0.1:1];
% % wmin = [0.1:0.1:1];
% tuning_fun_choose=[1, 3, 4, 5, 6, 7]; %function 2 was been deleted.
% tuning_runs = 5;

% row_counter = 1;

% % for w1 = wmax
% 	% for c1 = C1

% 		% if w1 > w2;
% 		% 	w1, w2
% 		% 	Inertia_Weight = w1 - (1:iter_max) .* (w1-w2) ./ iter_max;
% 			row_counter
% 			for i = 1:length(tuning_fun_choose)
% 				func_num = tuning_fun_choose(i);
% 				for j = 1:tuning_runs

%                     [gbest, gbestval, FES] =  PSO_func(fhd, D, pop_size, iter_max, Xmin, Xmax, Inertia_Weight, C_set, func_num);
					
% 					[gbest, gbestval, FES] =  CLPSO_func(fhd, D, pop_size, iter_max, Max_FES, Xmin, Xmax, Inertia_Weight, [c1, c1], func_num);
% 					xbest(j, :) = gbest;
% 	        		fbest(i, j) = gbestval;
% 				end
% 				f_mean(i) = mean(fbest(i, :));
% 				f_stdvar(i) = std(fbest(i, :));
% 			end
% 		% UGLY!!!!!!!!!!
% 			var_to_save(row_counter, 1) = c1;
% 			var_to_save(row_counter, 2) = c1;
% 			var_to_save(row_counter, 3) = f_mean(1);
% 			var_to_save(row_counter, 4) = f_stdvar(1);
% 			var_to_save(row_counter, 5) = f_mean(2);
% 			var_to_save(row_counter, 6) = f_stdvar(2);
% 			var_to_save(row_counter, 7) = f_mean(3);
% 			var_to_save(row_counter, 8) = f_stdvar(3);
% 			var_to_save(row_counter, 9) = f_mean(4);
% 			var_to_save(row_counter, 10) = f_stdvar(4);
% 			var_to_save(row_counter, 11) = f_mean(5);
% 			var_to_save(row_counter, 12) = f_stdvar(5);
% 			var_to_save(row_counter, 13) = f_mean(6);
% 			var_to_save(row_counter, 14) = f_stdvar(6);
% 			% for k = 1:6	
% 			% 	h = k + 1 
% 			% 	var_to_save(row_counter, k + 2) = f_mean(k);
% 			% 	var_to_save(row_counter, h + 2) = f_stdvar(k);
% 			% end
% 			row_counter = row_counter + 1;	
% 		% end
% 	end
% % end
% save('CLPSO_tuning_wmax_wmin.mat', 'var_to_save')


% % for i = 1:29
% %     func_num = i;
% %     for j = 1:runs
% %         i, j,
% %         [gbest, gbestval, FES] =  PSO_func(fhd, D, pop_size, iter_max, Xmin, Xmax, Inertia_Weight, C_set, func_num);
% %         xbest(j, :) = gbest;
% %         fbest(i, j) = gbestval;
% %         % fbest(i, j)
% %     end
% %     f_mean(i) = mean(fbest(i, :));
% % end
% % CLPSO_func(fhd, Dimension, Particle_Number, Max_Gen, Max_FES, VRmin, VRmax, Inertia_Weight, C_set, func_num)


% % for i=1:29
% % eval(['load input_data/shift_data_' num2str(i) '.txt']);
% % eval(['O=shift_data_' num2str(i) '(1:10);']);
% % f(i)=cec14_func(O',i);i,f(i)
% % end