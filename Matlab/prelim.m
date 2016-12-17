[ ipc,lifetime,energy,app_name ] = loaddata();
%%
Z{1} = ipc;
Z{2} = lifetime;
Z{3} = energy;
Y_name{1} = 'IPC';
Y_name{2} = 'Lifetime';
Y_name{3} = 'Energy';
% Parameters
n = size(Z{1},1);
m = length(app_name);
accuracy = zeros(m,2);

%% Sample
numSamples = 30;
id1 = 1:ceil(n/numSamples):n; % points uniform over 1:1024 
id2 = randperm(n); id2 = id2(1:numSamples); %random points 

%% Run
for Y_nameId = 1:length(Z),
    for i = 1:m
        ZZ = Z;
        [ acc, w_pred,online,offline ] = splitEM( [],ZZ,Y_nameId,id1,i );
        accuracy(i,Y_nameId) = acc;
        wl{Y_nameId,i} = w_pred;
        accuracy_offline(i,Y_nameId) = accuracy_rss(ZZ{Y_nameId}(:,i), offline);
        wl_offline{Y_nameId,i} = offline;
    end
end

mean(max(accuracy,0))
mean(max(accuracy_offline,0))
%% Plot
pp=3;
for Y_nameId = 1:length(Z),
    figure(Y_nameId);
    for i = 1:m    
        subplot(pp,ceil(m/pp),i);
        hold on;
        aa = Z{Y_nameId}(:,i); ind = 1:n;
        bb = cell2mat(wl{Y_nameId,i});
        %dd = wl_offline{Y_nameId,i};
        plot(aa); 
        plot(bb(ind));
        %plot(dd(ind));
        title([Y_name{Y_nameId},' - ',app_name{i},' : ',sprintf('%5.2f',accuracy(i,Y_nameId))]);
        xlim([1 n]);
    end
    %legend('True','EM','Offline');
    legend('True','EM');
end


%% Print config file
true_conf.cores = X(:,1); 
true_conf.freq  = X(:,2);
est_conf.cores = X(:,1); 
est_conf.freq  = X(:,2);
online_conf.cores = X(:,1); 
online_conf.freq  = X(:,2);
offline_conf.cores = X(:,1); 
offline_conf.freq  = X(:,2);
for i = 1:m      
    true_conf.perf  = Z{1}(:,i); 
    true_conf.lifetime = Z{2}(:,i);
    T = struct2table(true_conf);
    name = strcat('Poet_Config/true_',app_name(i),'.txt');
    writetable(T,name{1},'Delimiter',' ');

    est_conf.perf  = cell2mat(wl{1,i}); 
    est_conf.lifetime = cell2mat(wl{2,i});
    T = struct2table(est_conf);
    name = strcat('Poet_Config/est_',app_name(i),'.txt');
    writetable(T,name{1},'Delimiter',' ');
    
    online_conf.perf  = wl_online{1,i}; 
    online_conf.lifetime = wl_online{2,i};
    T = struct2table(online_conf);
    name = strcat('Poet_Config/online_',app_name(i),'.txt');
    writetable(T,name{1},'Delimiter',' ');
    
    offline_conf.perf  = wl_offline{1,i}; 
    offline_conf.lifetime = wl_offline{2,i};
    T = struct2table(offline_conf);
    name = strcat('Poet_Config/offline_',app_name(i),'.txt');
    writetable(T,name{1},'Delimiter',' ');
end

accuracy_out.name          = app_name;
accuracy_out.est_rate      = accuracy(:,1);
accuracy_out.est_rate      = accuracy(:,1);
accuracy_out.est_power     = accuracy(:,2);
accuracy_out.online_rate   = accuracy_online(:,1);
accuracy_out.online_power  = accuracy_online(:,2);
accuracy_out.offline_rate  = accuracy_offline(:,1);
accuracy_out.offline_power = accuracy_offline(:,2);
T = struct2table(accuracy_out);
name = strcat('Poet_Config/accuracy.txt');
writetable(T,name,'Delimiter',' ');
%% Nuclear norm
% addpath('MatrixCompletion/');
% lambda_tol = 10^(-6);
% tol        = 10^(-6);
% display    = 1;
% N          = 10;
% mode       = 'nuclear';
% A          = Z{1,1};
% for i = 1:m
% A(:,i) = A(:,i)/max(A(:,i));
% end
% 
% for i = 1:m
%     % Parameters
%     accuracy = zeros(m,2);
%     % Sample
%     numSamples = 20;
%     id1 = 1:ceil(n/numSamples):n; % points uniform over 1:1024 
%     id2 = randperm(n); id2 = id2(1:(n-numSamples)); %random points 
%     B = ones(size(A));
%     B(id2,i) = 0;
%     [CompletedMat, ier] = MatrixCompletion(A.*B, B,N, mode, lambda_tol, tol, display); 
%     wl_nuclear(:,i) = CompletedMat(:,i);
%     accuracy_norm(i) = 1 - norm(CompletedMat(:,i) - A(:,i))/norm(A(:,i));    
% end
% 
% figure(3);
% for i = 1:m    
%     subplot(pp,ceil(m/pp),i);
%     hold on;
%     plot(wl_nuclear(:,i));
%     plot(A(:,i),'r');    
%     title(['ipc - ',app_name{i},' : ',sprintf('%5.2f',accuracy_norm(i)) ]); 
%     hold off;
% end
