function [est_rate_pow,online_rate_pow,offline_rate_pow  ] = overhead()

[ rate,power,X,app_name ] = loaddata( 1 );
m = length(app_name);
Z{1} = rate;
Z{2} = power;
[n] = size(Z{1},1);
accuracy = zeros(m,2);

%% Sample
Sample_size = 1:128;

for j = 1:5
    
for k = 1:length(Sample_size)
    numSamples = Sample_size(k);
    %id1 = 1:ceil(n/numSamples):n; % points uniform over 1:1024 
    %id2 = randperm(n); id2 = id2(1:numSamples); %random points 
    id1 = randsample(n,numSamples);
    %% Run
    for Y_nameId = 1:2,
        for i = 1:m
            ZZ = Z;
            [ acc, w_pred,online,offline ] = splitEM( X,ZZ,Y_nameId,id1,i );
            accuracy(i,Y_nameId) = acc;
            wl{Y_nameId,i} = w_pred;
        
            accuracy_online(i,Y_nameId) = accuracy_rss(ZZ{Y_nameId}(:,i), online);
            accuracy_offline(i,Y_nameId) = accuracy_rss(ZZ{Y_nameId}(:,i), offline);
            wl_online{Y_nameId,i} = online;
            wl_offline{Y_nameId,i} = offline;
        end
    end
    est_rate_pow(j,k,:) = mean(max(accuracy,0));
    online_rate_pow(j,k,:) = mean(max(accuracy_online,0));
    offline_rate_pow(j,k,:) =  mean(max(accuracy_offline,0));
end
end
end

