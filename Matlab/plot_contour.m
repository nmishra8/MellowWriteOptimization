function [  ] = plot_contour( X,rate,titlee,ylab )
addpath('../MatlabCommon/');
addpath('../MatlabCommon/export_fig-master/');
coress_big = repmat([8,7,6,5,],19,1); 
cores_little = repmat([4,3,2,1],13,1);
X      = [[coress_big(:);cores_little(:)],X(:,2)];

x = unique(X(:,1));
y = unique(X(:,2));
z = zeros(length(x),length(y));
    for i = 1:length(x)
        for j = 1:length(y)
            index = find((X(:,1)==x(i)).*(X(:,2)==y(j)));
            if(~isempty(index));    z(i,j) = rate(index); end
        end
    end
    
    ax = contourf(x,y,z');
    %colormap(bone)
    colormap(brewermap([],'Blues'))


    title(titlee)
    ylabel('DVFS setting');
    %xlabel({'LITTLE                   big','Cores'});
    xlabel({'LITTLE Cores    big Cores'});
    
    %xlabel('Cores');
    
    %h = set(gca,'Xtick',x,'XTickLabel',{'0x1','0x3','0x7','0xF','0x10','0x30','0x70','0xF0',});
    fontsizee = 40;
    h = set(gca,'Xtick',x,'XTickLabel',{'1','2','3','4','1','2','3','4',},'FontName','Times New Roman');
    c = colorbar; % c.Label.String = ylab;c.Label.FontSize=fontsizee;
    set(gca,'fontsize',fontsizee);
    y1=get(gca,'ylim');
    
    cc=1;
%hold on
%plot([4.5 4.5],y1,'Color','k')

    %ax = gca;
    %ax.XTickLabelRotation = 45;
    
end

