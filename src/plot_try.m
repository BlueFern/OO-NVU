% %% Sort best results
% %sel_sum_prob(:,:) = zeros(length(best_results(:,1)),ini_NrFitVal);
% best_l = length(best_results(:,1));
%     for sel_j = 1 : ini_NrFitVal                                               %??
%         
%         for sel_i = 1:best_l                                
%             sel_prob(sel_i,sel_j)       = (abs(ini_FitVal(sel_j) - best_results(sel_i,sel_j+3)));     %probability for each individual to be chosen
%             if sel_prob(sel_i,sel_j) == 0;
%                 sel_prob(sel_i,sel_j) = 0.0000000000000000001;
%             end
%         end
%         %saving the probabilities
%         sel_sum(:,sel_j) = sum(sel_prob(:,sel_j));
%         
%         for sel_i = 1 : best_l
%             sel_prob_1(sel_i,sel_j) = sel_sum(sel_j) ./ sel_prob(sel_i,sel_j);
%         end
%         
%         sel_sum_2(:,sel_j) = sum(sel_prob_1(:,sel_j));
%         
%         for sel_i = 1 : best_l
%             sel_prob_2(sel_i,sel_j) = sel_prob_1(sel_i,sel_j) ./ sel_sum_2(:,sel_j);
%         end
%     end
%     
%     for ini_i = 1 : best_l
%         sel_prob_vec_1(ini_i,:) = sum(sel_prob_2(ini_i,:));
%         sel_prob_vec = sel_prob_vec_1 ./ ini_NrFitVal;
%     end
% clear sel_cum_probvec    
% sel_cum_probvec(:,:) = cumsum(sel_prob_vec);
% best_results  = [best_results(:,1), sel_cum_probvec, sel_prob_vec, best_results(:,4:sum(ini_NrPar)+3+ini_NrFitVal)];
% %% Plot Mutation rates 
%  
% figure(plot_k + 1)
% plot(linspace(1,r_simu_2-1,r_simu_2-1)',mu_p_total');
% legend('Change of average mutation rate');
% title('Change of average mutation rate');
% xlabel('# of generation');
% ylabel('probability of mutation');
% grid on;
% grid minor;
% 
% %% Plot best results
% plot_max = sortrows(best_results,3);
% result_string = {'Best result','Second best result','Third best result','4th best result','5th best result','6th best result','7th best result'};
% color_string = {'b','g','r','y','m','c','k'};
% plot_w = 1;
% plot_i = 1;
% plot_j = 1;
% plot_max_unique = zeros(7,length(best_results(1,:)));
% 
% while plot_w ==1 && plot_j <=7
%     if plot_max(best_l+1-plot_i,3) == plot_max(best_l-plot_i,3)
%         plot_w=1;
%     else
%        plot_max_unique(plot_j,:) = plot_max(best_l+1-plot_i,:);
%        plot_j = plot_j + 1;
%     end
%     plot_i = plot_i + 1;
% end
% 
% for plot_i = 1 : 7
% plot_max_line_1 = plot_max_unique(plot_i,:);
% x = 1:length(ini_ParamString);
% y_1 = plot_max_line_1(3+ini_NrFitVal+1 : 3+ini_NrFitVal+sum(ini_NrPar));
% plot(x, y_1, color_string{plot_i})
% hold on;
% end
% legend(result_string);
% set(gca, 'XTick',1:length(ini_ParamString), 'XTickLabel',ini_ParamString)
% grid on;
% grid minor;


time_beginning = clock;
%Start time of simulation
time_max = 10;
%Maximum simulation time in seconds
time_elapsed = 0;
%Just a preallocation

while time_elapsed < time_max
  
   
% HERE IS THE CODE

time_end     = clock;
time_elapsed = etime(time_end,time_beginning);
%Elapsed time between start of simulation and "now"
end

