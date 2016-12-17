%% phase change
clear; close all; 
addpath('../MatlabCommon/');
addpath('../MatlabCommon/export_fig-master/');

figure(1);
xx = load('../results/STREAM-phases/leo.txt');
subplot(2,1,1); hold on; plot(xx(:,2));title('STREAM');
subplot(2,1,2); hold on; plot(xx(:,3));

xx = load('../results/STREAM-phases/leopoet.txt');
subplot(2,1,1); hold on; plot(xx(:,2));
subplot(2,1,2); hold on; plot(xx(:,3));

xx = load('../results/STREAM-phases/poet.txt');
subplot(2,1,1); hold on; plot(xx(:,2));
subplot(2,1,2); hold on; plot(xx(:,3));
legend('leo','leopoet','poet');

figure(2);
xx = load('../results/x264-native-ducks/leo.txt');
subplot(2,1,1); hold on; plot(xx(:,2));title('X264')
subplot(2,1,2); hold on; plot(xx(:,3));

xx = load('../results/x264-native-ducks/leopoet.txt');
subplot(2,1,1); hold on; plot(xx(:,2));
subplot(2,1,2); hold on; plot(xx(:,3));

xx = load('../results/x264-native-ducks/poet.txt');
subplot(2,1,1); hold on; plot(xx(:,2));
subplot(2,1,2); hold on; plot(xx(:,3));
legend('leo','leopoet','poet');

% better results

legend_name = {'leo','poet','leo-poet'};
color_scheme = load('color_scheme.mat');
color_scheme = color_scheme.abcd;
color_scheme = color_scheme.col(ismember(color_scheme.name,legend_name),:);
        
subplot(2,1,1); 
hold on; 
plot(Frame,LEO_PERF);%,'Color',color_scheme(1,:));
plot(Frame,POET_PERF);
plot(Frame,LEOPOET_PERF);
hold off;
colormap(color_scheme);
title('Performance for X264');

subplot(2,1,2); 
hold on; 
plot(Frame,LEO_POWER);%,'Color',color_scheme(1,:));
plot(Frame,POET_POWER);
plot(Frame,LEOPOET_POWER);
hold off;
colormap(color_scheme)
legend(legend_name);
title('Power for X264');

%% single app

[ rate,power,X,app_name ] = loaddata(0 );
m  = length(app_name);
pp = 4;
Z{1} = rate;
Z{2} = power;
% scatter power vs performance
for i = 1:m    
    figure(3);
    subplot(pp,ceil(m/pp),i);
    hold on;    
    scatter(Z{2}(:,i),Z{1}(:,i)); 
    hold off;
    title([app_name{i}]); 
    ylabel('rate');
    xlabel('power');
end

% bar plot comparison
% legend_name = {'optimal','online','offline','leo','poet','leo-poet'};
% color_scheme = brewermap(length(legend_name)+1,'BrBG');   
% color_scheme(sum(color_scheme,2)==max(sum(color_scheme,2)),:)=[];
% abcd.name = legend_name';
% abcd.col = color_scheme;
% abcd = struct2table(abcd);
close all;
legend_name = {'optimal','online','offline','leo','poet','leo-poet'};
xx_enr = barplot_target( '../results/single-app/dyn-eff-', 'energy',legend_name);
%export_fig ~/Desktop/LEO-Bard/Matlab/Figure/single-app-energy.pdf -transparent -painters
xx_err = barplot_target( '../results/single-app/dyn-mape-', 'performance', legend_name);
%export_fig ~/Desktop/LEO-Bard/Matlab/Figure/single-app-mape.pdf -transparent -painters
close all;
legend_name = {'leo','poet','leo-poet'};
xx_enr_m = barplot_target( '../results/multi-app/ma-eff-', 'energy',legend_name);
%export_fig ~/Desktop/LEO-Bard/Matlab/Figure/multi-app-energy.pdf -transparent -painters
xx_err_m = barplot_target( '../results/multi-app/ma-err-', 'performance',legend_name);
%export_fig ~/Desktop/LEO-Bard/Matlab/Figure/multi-app-mape.pdf -transparent -painters

%% scatter
legend_name = {'optimal','online','offline','leo','poet','leo-poet'};
target_perf = [0.1,0.3,0.5,0.7,0.9];
n = length(target_perf);
color_scheme = load('color_scheme.mat');
color_scheme = color_scheme.abcd;
color_scheme = color_scheme.col(ismember(color_scheme.name,legend_name),:);
close all;
for i = 1: n
    filename = '../results/single-app/dyn-eff-';
    file_name = [filename,num2str(target_perf(i)),'.txt'];
    fid = fopen(file_name);
    c = textscan(fid,'%s %f %f %f %f %f %f','HeaderLines', 0 );
    fclose(fid);
    xx_err = max(cell2mat(c(2:end)),0);
    xx_err(:,find(any(~isnan(xx_err),1)==0))=[];
    
    filename = '../results/single-app/dyn-mape-';
    file_name = [filename,num2str(target_perf(i)),'.txt'];
    fid = fopen(file_name);
    c = textscan(fid,'%s %f %f %f %f %f %f','HeaderLines', 0 );
    fclose(fid);
    app_name = c{1};
    xx_enr = max(cell2mat(c(2:end)),0);
    xx_enr(:,find(any(~isnan(xx_err),1)==0))=[];
    
    subplot(n,1,i); 
    for j = 1:length(legend_name)
        hold on;
        scatter(xx_err(:,j),xx_enr(:,j),'MarkerEdgeColor',color_scheme(j,:),...
              'MarkerFaceColor',color_scheme(j,:),'LineWidth',1.5);
        hold off;
    end
   
    if(i==1)
        l = legend(legend_name,'Location','north','Orientation','horizontal');
        set(l,'interpreter', 'none');
        %title(titlee);
    end
%     if(i==n)
%         xlabel('APP NAME');
%     else
%         set(gca, 'XTick', []);
%     end
%     
%     ylim([min(min(xx_err)),max(max(xx_err))]);
%     xlim([0 (length(app_name)+1)]);
    %set(gca,'fontsize',13);
    
   % ylabel(num2str(target_perf(i)));
     ylabel('MAPE');
    xlabel('Normalized energy (higher is better)');
    grid on;
end

legend_name = {'leo','poet','leo-poet'};
%legend_name = {'leo','poet','leo-poet'};
target_perf = [0.1,0.3,0.5,0.7,0.9];
n = length(target_perf);
color_scheme = load('color_scheme.mat');color_scheme = color_scheme.abcd;
color_scheme = color_scheme.col(ismember(color_scheme.name,legend_name),:);
close all;
ind = 1;
for i = 1: n
    filename = '../results/multi-app/ma-eff-';
    file_name = [filename,num2str(target_perf(i)),'.txt'];
    fid = fopen(file_name);
    c = textscan(fid,'%s %f %f %f','HeaderLines', 0 );
    fclose(fid);
    xx_err = max(cell2mat(c(2:end)),0);
    xx_err(:,find(any(~isnan(xx_err),1)==0))=[];
    
    filename = '../results/multi-app/ma-err-';
    file_name = [filename,num2str(target_perf(i)),'.txt'];
    fid = fopen(file_name);
    c = textscan(fid,'%s %f %f %f','HeaderLines', 0 );
    fclose(fid);
    app_name = c{1};
    xx_enr = max(cell2mat(c(2:end)),0);
    xx_enr(:,find(any(~isnan(xx_err),1)==0))=[];
    
     
    for j = 1:length(legend_name)
        %hold on;
        
        subplot(n,3,ind);
        scatter(xx_err(:,j),xx_enr(:,j),'MarkerEdgeColor',color_scheme(j,:),...
              'MarkerFaceColor',color_scheme(j,:),'LineWidth',1.5);
        ind = ind+1;
        %hold off;
    end
   
    if(i==1)
        l = legend(legend_name,'Location','north','Orientation','horizontal');
        set(l,'interpreter', 'none');
        title(['multi-app ',num2str(target_perf(i))]);
    end
%     if(i==n)
%         xlabel('APP NAME');
%     else
%         set(gca, 'XTick', []);
%     end
%     
%     ylim([min(min(xx_err)),max(max(xx_err))]);
%     xlim([0 (length(app_name)+1)]);
    %set(gca,'fontsize',13);
    
    ylabel('MAPE');
    xlabel('Normalized energy (higher is better)');
    grid on;
end

%%  heatmap for all
clear;close all;
[ rate,power,X,app_name ] = loaddata( 1 );
for k = 1:size(rate,2)
    subplot(5,5,k);
    plot_contour( X,rate(:,k),app_name(k));
end

close all;
[ rate,power,X,app_name ] = loaddata( 1 );
for k = 1:size(rate,2)
    subplot(5,5,k);
    plot_contour( X,power(:,k),app_name(k));
end
% heatmap for x264

str = 'X264_ducks';
filename = strcat('~/Dropbox/Poet_Config/true_',str,'.txt');
true = importdata(filename, ' ', 1);true = true.data;
filename = strcat('~/Dropbox/Poet_Config/est_',str,'.txt');
est = importdata(filename, ' ', 1); est = est.data;
filename = strcat('~/Dropbox/Poet_Config/online_',str,'.txt');
online = importdata(filename, ' ', 1); online = online.data;
filename = strcat('~/Dropbox/Poet_Config/offline_',str,'.txt');
offline = importdata(filename, ' ', 1); offline = offline.data;

close all;
subplot(2,4,1); plot_contour( X,true(:,2),[str,' True', 'Performance']);
subplot(2,4,2); plot_contour( X,est(:,2),[str,' LEO', 'Performance']);
subplot(2,4,3); plot_contour( X,online(:,2),[str,' Online', 'Performance']);
subplot(2,4,4); plot_contour( X,offline(:,2),[str,' Offline', 'Performance']);

subplot(2,4,5); plot_contour( X,true(:,3),[str,' True', 'Power']);
subplot(2,4,6); plot_contour( X,est(:,3),[str,' LEO', 'Power']);
subplot(2,4,7); plot_contour( X,online(:,3),[str,' Online', 'Power']);
subplot(2,4,8); plot_contour( X,offline(:,3),[str,' Offline', 'Power']);

% heatmap for KMEANS LUD X264
close all;
subplot(2,3,1);plot_contour( X,rate(:,12),app_name(12),'Performance (in Heartbeat/sec)');
subplot(2,3,2);plot_contour( X,rate(:,14),app_name(14),'Performance (in Heartbeat/sec)');
subplot(2,3,3);plot_contour( X,rate(:,23),'X264-native','Performance (in Heartbeat/sec)');
subplot(2,3,4);plot_contour( X,power(:,19),app_name(19),'Power (in Watts)');
subplot(2,3,5);plot_contour( X,power(:,16),app_name(16),'Power (in Watts)');
subplot(2,3,6);plot_contour( X,power(:,23),'X264-native','Power (in Watts)');

%export_fig ~/Desktop/ASPLOS-17/figures/sample-contour3.png -transparent -painters
%export_fig ~/Desktop/LEO-Bard/Matlab/Figure/performance-contour3.pdf -transparent -painters
%export_fig ~/Desktop/ASPLOS-17/figures/performance-contour3.pdf
%-transparent -painters -opengl
%subplot(2,3,4);plot_contour( X,power(:,12),app_name(12));
%subplot(2,3,5);plot_contour( X,power(:,16),app_name(16));
%subplot(2,3,6);plot_contour( X,power(:,23),'X264-native');


%%
plot(est_rate_pow(:,1))
hold on;
plot(offline_rate_pow(:,1))

%% plot for memory and clockspeed
clear;close all;
[ rate,power,X,app_name ] = loaddata( 0 );
rate_4 = rate(1:19,:); 
x = X(1:19,2);
R2summary=[];
for k = 1:size(rate,2)
    subplot(5,5,k);
    scatter(x,rate_4(:,k));
    regs = regstats(x,rate_4(:,k),'linear');
    R2summary(end+1) = regs.adjrsquare; 
    xlabel('Clockspeed');
    ylabel('Rate');
    title(app_name(k));
end
