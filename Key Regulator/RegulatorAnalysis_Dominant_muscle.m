%% Load files
load('./Resource/Regulator-DominantMetabolite_muscle.mat');

%% Parameters - Metabolite
temp_Info = Info_Metabolite;
temp_data = data_Metabolite;
temp_xlabel = 'Differential Metabolites';
temp_ylabel = 'Number of Regulations';

name_figure = 'RegulatorAnalysis-DominantMetabolite';
bar_YLim = [0,300];
bar_FontSize = 8;
cutoff = 30;

%% Parameters - Metabolite_SP
temp_Info = Info_Metabolite_SP;
temp_data = data_Metabolite_SP;
temp_xlabel = 'Differential Metabolites (Substrate & Product)';
temp_ylabel = 'Number of Regulations';

name_figure = 'RegulatorAnalysis-DominantMetabolite_SP';
bar_YLim = [0,100];
bar_FontSize = 8;
cutoff = 10;

%% Parameters - Metabolite_CF
temp_Info = Info_Metabolite_CF;
temp_data = data_Metabolite_CF;
temp_xlabel = 'Differential Metabolites (Cofactor)';
temp_ylabel = 'Number of Regulations';

name_figure = 'RegulatorAnalysis-DominantMetabolite_CF';
bar_YLim = [0,100];
bar_FontSize = 8;
cutoff = 1;

%% Parameters - Metabolite_AC
temp_Info = Info_Metabolite_AC;
temp_data = data_Metabolite_AC;
temp_xlabel = 'Differential Metabolites (Activator)';
temp_ylabel = 'Number of Regulations';

name_figure = 'RegulatorAnalysis-DominantMetabolite_AC';
bar_YLim = [0,15];
bar_FontSize = 8;
cutoff = 6;

%% Parameters - Metabolite_IN
temp_Info = Info_Metabolite_IN;
temp_data = data_Metabolite_IN;
temp_xlabel = 'Differential Metabolites (Inhibitor)';
temp_ylabel = 'Number of Regulations';

name_figure = 'RegulatorAnalysis-DominantMetabolite_IN';
bar_YLim = [0,100];
bar_FontSize = 8;
%cutoff = 10;
cutoff = 17;

%% Plot
temp_include = temp_data(:,1) >= cutoff;
temp_Info = temp_Info(temp_include,:);
temp_data = temp_data(temp_include,:);

bar_data = temp_data(:,1);
bar_XTick = 1:1:size(temp_Info,1);
bar_XTickLabel = cell(size(temp_Info,1),1);
for i = 1:size(temp_Info,1)
    temp_name = temp_Info{i,2};
    temp_color = temp_Info{i,3};
    if strcmp(temp_color,'#CCCCFF')
        temp_prefix = '\color[rgb]{0.6,0.6,1}';
    elseif strcmp(temp_color,'#FFCCCC')
        temp_prefix = '\color[rgb]{1,0.6,0.6}';
    elseif strcmp(temp_color,'#FFCCFF')
        temp_prefix = '\color[rgb]{1,0.6,1}';
    end
    bar_XTickLabel{i} = [temp_prefix,temp_name];
end

f = figure(1);
hold on;
index_list = 1:sum(temp_include);

temp_index_list_WT = strcmp(temp_Info(temp_include,3),'#CCCCFF');
b1 = bar(index_list(temp_index_list_WT),bar_data(temp_index_list_WT),'FaceColor','#CCCCFF');
temp_index_list_OB = strcmp(temp_Info(temp_include,3),'#FFCCCC');
b2 = bar(index_list(temp_index_list_OB),bar_data(temp_index_list_OB),'FaceColor','#FFCCCC');
temp_index_list_MIX = strcmp(temp_Info(temp_include,3),'#FFCCFF');
b3 = bar(index_list(temp_index_list_MIX),bar_data(temp_index_list_MIX),0.8,'FaceColor','#FFCCFF');

xtips1 = b1.XEndPoints;
ytips1 = b1.YEndPoints;
labels1 = string(b1.YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',bar_FontSize);
xtips2 = b2.XEndPoints;
ytips2 = b2.YEndPoints;
labels2 = string(b2.YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',bar_FontSize);
xtips3 = b3.XEndPoints;
ytips3 = b3.YEndPoints;
labels3 = string(b3.YData);
text(xtips3,ytips3,labels3,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',bar_FontSize);

xlabel(temp_xlabel);
ylabel(temp_ylabel);

ax = gca(f);
ax.XTick = bar_XTick;
ax.XTickLabel = bar_XTickLabel;
ax.XTickLabelRotation = -45;  %%%%%
ax.YLim = bar_YLim;
legend("WT-higher","\it{ob}\rm{/}\it{ob}\rm{-higher}","mixed pattern");

%% Save
exportgraphics(f,['./Fig/',name_figure,'_muscle.pdf'],'ContentType','vector');

%%
