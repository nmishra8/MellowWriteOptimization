function [ ] = linearProgramm()


phases = 0 ;
if(phases == 0)
    %app_name = {'BACKPROP','BFS','BLACKSCHOLES_10M','BODYTRACK','FACESIM','FERRET',...
     %       'HEARTWALL', 'HOTSPOT','JACOBI','KMEANS','KMEANSNF','LAVAMD',...
      %      'LEUKOCYTE','LUD','NW','SHA',...
       %     'SRAD','STREAM','STREAM_threads','X264_ducks','X264_native'}';
       app_name = {'KMEANS','KMEANSNF','LAVAMD','LEUKOCYTE'}';
else
    app_name = {'BACKPROP','BFS','BLACKSCHOLES_1M','BLACKSCHOLES_10M','BLACKSCHOLES_10M_1','BODYTRACK','FACESIM','FERRET',...
            'HEARTWALL', 'HOTSPOT','JACOBI','KMEANS','KMEANSNF','LAVAMD',...
            'LEUKOCYTE','LUD','NW','SHA',...
            'SRAD','STREAM','STREAM_threads','X264_ducks','X264_native','X264_phases'}';
end  



for i= 1:length(app_name)
str = app_name(i);
filename = strcat('~/Dropbox/Poet_Config/true_',str,'.txt');
true = importdata(filename{1}, ' ', 1);true = true.data;
filename = strcat('~/Dropbox/Poet_Config/est_',str,'.txt');
est = importdata(filename{1}, ' ', 1); est = est.data;
filename = strcat('~/Dropbox/Poet_Config/online_',str,'.txt');
online = importdata(filename{1}, ' ', 1); online = online.data;
filename = strcat('~/Dropbox/Poet_Config/offline_',str,'.txt');
offline = importdata(filename{1}, ' ', 1); offline = offline.data;


w = 0.1:0.1:max(true(:,2));
trueEnergy = zeros(size(w));
estEnergy = zeros(size(w));
onlineEnergy = zeros(size(w));
offlineEnergy = zeros(size(w));

T = 1.5;
Aeq = ones(1, 128);
beq = T;
lb = zeros(128 ,1);
ub = T*ones(128 ,1);
options=optimset('Display', 'off');
for j = 1:length(w)
% For true data
        f = true(:,3);
        A = -true(:,2)';
        b = -w(j);
        [~,fval,~,~,~] = linprog(f,A,b,Aeq,beq,lb,ub,[],options);
        trueEnergy(j) = fval; 
        
        % LEO
        
        f = est(:,3);
        A = [-true(:,2)';-est(:,2)'];
        b = [-w(j),-w(j)];
        [xest,fval,~,~,~] = linprog(f,A,b,Aeq,beq,lb,ub,[],options);
        estEnergy(j) = true(:,3)'*xest; 
        
        % online
        f = online(:,3);
        A = [-true(:,2)';-online(:,2)'];
        b = [-w(j),-w(j)];
        [xest,fval,~,~,~] = linprog(f,A,b,Aeq,beq,lb,ub,[],options);
        onlineEnergy(j) = true(:,3)'*xest; 
        
        % offline
        f = offline(:,3);
        A = [-true(:,2)';-offline(:,2)'];
        b = [-w(j),-w(j)];
        [xest,fval,~,~,~] = linprog(f,A,b,Aeq,beq,lb,ub,[],options);
        offlineEnergy(j) = true(:,3)'*xest; 
        
end
fprintf('Perf acc - EST: %f, ONLINE: %f, OFFLINE: %f \n',...
    accuracy_rss(true(:,2),est(:,2)), accuracy_rss(true(:,2),online(:,2)),...
    accuracy_rss(true(:,2),offline(:,2)));
fprintf('Power acc - EST: %f, ONLINE: %f, OFFLINE: %f\n',...
    accuracy_rss(true(:,3),est(:,3)), accuracy_rss(true(:,3),online(:,3)),...
    accuracy_rss(true(:,3),offline(:,3)));

figure
subplot(1,3,1); 
hold on; 
plot(true(:,2),'linewidth',1.5);
plot(est(:,2),'linewidth',1.5);
plot(online(:,2),'linewidth',1.5);
plot(offline(:,2),'linewidth',1.5);
hold off;
set(gca,'fontsize',13);
title(str);
xlabel('Configuration Index');
ylabel('Performance rate (in iter/sec)');
subplot(1,3,2); 
hold on; 
plot(true(:,3),'linewidth',1.5);
plot(est(:,3),'linewidth',1.5);
plot(online(:,3),'linewidth',1.5);
plot(offline(:,3),'linewidth',1.5);
hold off;
set(gca,'fontsize',13);
title(str);
xlabel('Configuration Index');
ylabel('System-Power (in Watts)');
subplot(1,3,3); 
hold on; 
plot(w,trueEnergy,'linewidth',1.5);
plot(w,estEnergy,'linewidth',1.5);
plot(w,onlineEnergy,'linewidth',1.5);
plot(w,offlineEnergy,'linewidth',1.5);
hold off;
set(gca,'fontsize',13);
title(str);
xlabel('Utilization');
ylabel('Energy in (J)');
legend('True','LEO','Online','Offline','Location','northoutside','Orientation','horizontal');
%export_fig ~/Desktop/LEO-Bard/Matlab/Figure/X264_phases_performance.pdf -transparent -painters
%close all;

end

end
