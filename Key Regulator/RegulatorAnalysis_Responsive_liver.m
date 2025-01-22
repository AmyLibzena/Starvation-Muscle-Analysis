%% Load files
load('./Resource/Regulator-ResponsiveMetabolite_liver.mat');

%% Parameters - Metabolite
temp_Info = Info_Metabolite;
temp_data = data_Metabolite;
temp_xlabel = 'Responsive Metabolites';
temp_ylabel = 'Number of Regulations';

name_figure = 'RegulatorAnalysis-ResponsiveMetabolite';
bar_YLim = [-500,500];
bar_FontSize = 8;
cutoff = 30;

%% Plot
temp_include = temp_data(:,3) >= cutoff;
temp_Info = temp_Info(temp_include,:);
temp_data = temp_data(temp_include,:);

bar_data_WT = temp_data(:,1);
bar_data_OB = temp_data(:,2)*(-1);
bar_XTick = 1:1:size(temp_Info,1);
bar_XTickLabel = cell(size(temp_Info,1),1);
for i = 1:size(temp_Info,1)
    temp_name = temp_Info{i,2};
    temp_color = temp_Info{i,3};
    if strcmp(temp_color,'#00FF00')
        temp_prefix = '\color[rgb]{0,1,0}';
    elseif strcmp(temp_color,'#0000FF')
        temp_prefix = '\color[rgb]{0,0,1}';
    elseif strcmp(temp_color,'#FF0000')
        temp_prefix = '\color[rgb]{1,0,0}';
    elseif strcmp(temp_color,'#FF00FF')
        temp_prefix = '\color[rgb]{1,0,1}';
    end
    bar_XTickLabel{i} = [temp_prefix,temp_name];
end

f = figure(1);
hold on;
b1 = bar(bar_data_WT,'b');
b2 = bar(bar_data_OB,'r');

xtips1 = b1.XEndPoints;
ytips1 = b1.YEndPoints;
labels1 = string(b1.YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',bar_FontSize);

xtips1 = b2.XEndPoints;
ytips1 = b2.YEndPoints;
labels1 = string(b2.YData*(-1));
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','top','FontSize',bar_FontSize);

xlabel(temp_xlabel);
ylabel(temp_ylabel);

ax = gca(f);
ax.XTick = bar_XTick;
ax.XTickLabel = bar_XTickLabel;
ax.XTickLabelRotation = -45;  %%%%%
ax.YLim = bar_YLim;
legend("WT","\it{ob}\rm{/}\it{ob} ");

%% Save
exportgraphics(f,['./Fig/',name_figure,'_liver.pdf'],'ContentType','vector');

%%
