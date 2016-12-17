function [ xx_err ] = barplot_target( filename, titlee, legend_name)
% bar plot comparison
target_perf = [0.1,0.3,0.5,0.7,0.9];
n = length(target_perf);

color_scheme = load('color_scheme.mat');
color_scheme = color_scheme.abcd;
color_scheme = color_scheme.col(ismember(color_scheme.name,legend_name),:);

% energy
for i = 1: n
    file_name = [filename,num2str(target_perf(i)),'.txt'];
    fid = fopen(file_name);
    c = textscan(fid,'%s %f %f %f %f %f %f','HeaderLines', 0 );
    fclose(fid);
    
    app_name = c{1};
    xx_err = max(cell2mat(c(2:end)),0);
    xx_err(:,find(any(~isnan(xx_err),1)==0))=[];
    
    subplot(n,1,i); 
    bar(xx_err,'LineWidth',0.01,'EdgeColor',[1 1 1]);
    %colormap(color_scheme);
    if(i==1)
        l = legend(legend_name,'Location','north','Orientation','horizontal');
        set(l,'interpreter', 'none');
        title(titlee);
    end
    if(i==n)
        set(gca,'XTickLabel',app_name,'XTick',1:numel(app_name));
        xlabel('APP NAME');
    else
        set(gca, 'XTick', []);
    end
    
    ylim([min(min(xx_err)),max(max(xx_err))]);
    xlim([0 (length(app_name)+1)]);
    %set(gca,'fontsize',13);
    
    ylabel(num2str(target_perf(i)));
    
    grid on;
    ax = gca;
    ax.XTickLabelRotation = 45; 
end



end

