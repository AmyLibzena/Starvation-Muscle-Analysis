%% Load files
load('./Resource/Regulator-Metabolite.mat');
load('./Resource/Regulator-TF_ANOVA.mat');

%% Parameters - Metabolite
temp_Info = Info_Metabolite;
temp_data = data_Metabolite;
temp_xlabel = 'Metabolite';
temp_ylabel = 'Number of Regulations';

name_figure = 'RegulatorAnalysis-Metabolite';
bar_YLim = [-500,500];
bar_FontSize = 8;
cutoff = 20;

%% Parameters - Metabolite_SP
temp_Info = Info_Metabolite_SP;
temp_data = data_Metabolite_SP;
temp_xlabel = 'Metabolite (Substrate & Product)';
temp_ylabel = 'Number of Regulations';

name_figure = 'RegulatorAnalysis-Metabolite_SP';
bar_YLim = [-150,150];
bar_FontSize = 8;
cutoff = 10;

%% Parameters - Metabolite_CF
temp_Info = Info_Metabolite_CF;
temp_data = data_Metabolite_CF;
temp_xlabel = 'Metabolite (Cofactor)';
temp_ylabel = 'Number of Regulations';

name_figure = 'RegulatorAnalysis-Metabolite_CF';
bar_YLim = [-150,150];
bar_FontSize = 8;
cutoff = 1;

%% Parameters - Metabolite_AC
temp_Info = Info_Metabolite_AC;
temp_data = data_Metabolite_AC;
temp_xlabel = 'Metabolite (Activator)';
temp_ylabel = 'Number of Regulations';

name_figure = 'RegulatorAnalysis-Metabolite_AC';
bar_YLim = [-80,80];
bar_FontSize = 8;
cutoff = 5;

%% Parameters - Metabolite_IN
temp_Info = Info_Metabolite_IN;
temp_data = data_Metabolite_IN;
temp_xlabel = 'Metabolite (Inhibitor)';
temp_ylabel = 'Number of Regulations';

name_figure = 'RegulatorAnalysis-Metabolite_IN';
bar_YLim = [-150,150];
bar_FontSize = 8;
%cutoff = 10;
cutoff = 20;

%% Parameters - TF2Enzyme_T
temp_Info = Info_TF2Enzyme_T;
temp_data = data_TF2Enzyme_T;
temp_xlabel = 'Responsive TFs';
temp_ylabel = 'Number of Regulated Enzymes';

name_figure = 'RegulatorAnalysis-TF2Enzyme_T_ANOVA';
bar_YLim = [-300,300];
bar_FontSize = 8;
cutoff = 1;

%% Parameters - TF2Transporter_T
temp_Info = Info_TF2Transporter_T;
temp_data = data_TF2Transporter_T;
temp_xlabel = 'Responsive TFs';
temp_ylabel = 'Number of Regulated Transporters';

name_figure = 'RegulatorAnalysis-TF2Transporter_T_ANOVA';
bar_YLim = [-80,80];
bar_FontSize = 8;
cutoff = 1;

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
exportgraphics(f,['./Fig/',name_figure,'.pdf'],'ContentType','vector');

%%
