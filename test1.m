FootDropModel.simulate(1, 315, 12);

% generate curve fitted activation plots
% x = 0:0.01:1;
% resultS = zeros(length(x));
% resultTA = zeros(length(x));
% 
% for i = 1:length(x)
%     resultS(i) = HillTypeMuscle.SActivationRegression.eval(x(i));
%     resultTA(i) = HillTypeMuscle.TAActivationRegression.eval(x(i));
% end
% 
% figure
% hold on
% plot(x, resultS, 'r')
% title('Soleus Activation')
% ylabel('Activation')
% xlabel('Gait Cycle Completion')
% hold off
% figure
% hold on
% plot(x, resultTA, 'g')
% title('Tibialis Activation')
% ylabel('Activation')
% xlabel('Gait Cycle Completion')

