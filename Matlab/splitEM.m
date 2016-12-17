function [ accuracy, w_pred,online, offline] = splitEM( X,Z,Y_nameId,id1,i  )
Y_name{1} = 'IPC';
Y_name{2} = 'Lifetime';
Y_name{3} = 'Energy';
[n,m] = size(Z{1});
W = zeros(n,1);
W(id1) = 1;
numSamples = sum(W)
model = 'Quadratic';
whichstats = {'yhat','r','mse','rsquare','adjrsquare','beta','leverage','cookd','tstat'};
pi = 1;
tau = 1;
epsilon = 50000;
iteration_limit = 20;

    u = Z{Y_nameId}; 
	for j = 1:m;    y{j} = u(:,j);  end

    %% initializing the vector with missing values
    Y_known = W.*y{i};
    y_em = y; 
    y_em{i} = Y_known;
    I = diag(W.*ones(n,1));
    %% Start EM algorithm 
    fprintf('i = %f\n',i);
    SupportData = cell2mat(y); 
    SupportData(:,i) = [];
    
    % Offline initialization
    Old.mu = mean(SupportData,2);    
    offline = Old.mu;
    
    % Online initialization
    try
        Y_sample = y{i}(id1); 
        X_sample = X(id1,:);
        stats = regstats(Y_sample,X_sample,model,whichstats);
        online = x2fx(X, model)*stats.beta; 
    catch 
        online = zeros(size(offline));
    end
    
    %
    Old.C = speye(n);
    Old.sigma = 1; 
    startMain = tic;   
    [return_em, it] = runEM4(Old, W, cell2mat(y_em), i ,epsilon, iteration_limit, y);
     elapsedTime = toc(startMain);
    w_pred{i} = return_em.w;
    residual = w_pred{i} - y{i};
    rss = sum(residual.*residual); 
    tss = sum((y{i}-mean(y{i})).*(y{i}-mean(y{i})));
    residualsquare = 1-rss/tss;
    accuracy = 1-(1-residualsquare)*(n-1)/(n-2);
    
%     %% Plot
%     close all;
%     fontsize = 20;
%     h = figure;
%     hold on;
%     plot(y{i},'r','LineWidth',2.5);
%     plot(w_pred{i},'b','LineWidth',2.5); 
%     scatter(1:length(Y_known),Y_known,'c','fill');
%     hl = legend('Data','Estimate');
%     set(hl, 'FontSize',fontsize);
%     ylabel( Y_name{Y_nameId},'FontSize', fontsize);
%     xlabel('Configuration index','FontSize', fontsize);
%     annotation('textbox',[0.15 0.8 0.1 0.1], 'String', strcat('Accuracy = ',... 
%                 sprintf('%5.2f\n',accuracy)), 'FontSize', fontsize);
%     hold off;
    
   
    fprintf('accuracy = %f, time = %f\n\n', accuracy, elapsedTime);


end

