%% Load files
load('./Resource/ProteomeData_liver.mat');

%% Parameters
num_Genotype = 2;
index_WT = 1;
index_OB = 2;
num_TissueType = 1;
index_liver = 1;
num_TimePoint = 8;
TimePoint = [0,2,4,6,8,12,16,24];
col_num_liver = [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5];
num_Replicate = 5;
Group = {'WT0','WT0','WT0','WT0','WT0','WT2','WT2','WT2','WT2','WT2','WT4','WT4','WT4','WT4','WT4','WT6','WT6','WT6','WT6','WT6','WT8','WT8','WT8','WT8','WT8','WT12','WT12','WT12','WT12','WT12','WT16','WT16','WT16','WT16','WT16','WT24','WT24','WT24','WT24','WT24','OB0','OB0','OB0','OB0','OB0','OB2','OB2','OB2','OB2','OB2','OB4','OB4','OB4','OB4','OB4','OB6','OB6','OB6','OB6','OB6','OB8','OB8','OB8','OB8','OB8','OB12','OB12','OB12','OB12','OB12','OB16','OB16','OB16','OB16','OB16','OB24','OB24','OB24','OB24','OB24'};
threshold_q = 0.1;

%%% Check Enough Data (Exist in >= half replicates)
is_EnoughData = false(num_TissueType,size(Info_ProteomeName,1));
for i = 1:size(is_EnoughData,2)
    logic_EnoughData = CheckEnoughData(Proteome_liver(:,i), col_num_liver);
    is_EnoughData(index_liver,i) = logic_EnoughData;
end

%%% Calculation
%%% Mean
mean_liver = NaN(length(col_num_liver),size(Info_ProteomeName,1));
for i = 1:size(mean_liver,2)
    mean_liver(:,i) = GetMean(Proteome_liver(:,i),col_num_liver);
end

%%% FC & log2FC (WT tX/t0, OB tX/t0, OB/WT tX)(x = 0,2,4,6,8,12,16,24)
FC_liver = [mean_liver(1:num_TimePoint,:)./mean_liver(1,:);mean_liver(num_TimePoint+1:num_TimePoint*2,:)./mean_liver(num_TimePoint+1,:);mean_liver(num_TimePoint+1:num_TimePoint*2,:)./mean_liver(1:num_TimePoint,:)];
log2FC_liver = log2(FC_liver);
 
%%% p values (One-way ANOVA)(WT t0-tX;OB t0-tX;OB/WT tX)
p_liver = NaN(num_Genotype,size(Info_ProteomeName,1));
temp_num = sum(col_num_liver(1:num_TimePoint));
% liver - WT
for i = 1:size(Info_ProteomeName,1)
    p_liver(index_WT,i) = anova1(Proteome_liver(1:temp_num,i),Group(1:temp_num),'off');
end
% liver - OB
for i = 1:size(Info_ProteomeName,1)
    p_liver(index_OB,i) = anova1(Proteome_liver(temp_num+1:end,i),Group(temp_num+1:end),'off');
end
% liver - OB/WT
for i = 1:size(Info_ProteomeName,1)
    for j = 1:num_TimePoint
        if j == 1
            temp_list_WT = Proteome_liver(1:col_num_liver(j),i);
            temp_list_OB = Proteome_liver(temp_num+1:temp_num+col_num_liver(num_TimePoint+j),i);
        else
            temp_list_WT = Proteome_liver(sum(col_num_liver(1:j-1))+1:sum(col_num_liver(1:j)),i);
            temp_list_OB = Proteome_liver(sum(col_num_liver(1:num_TimePoint+j-1))+1:sum(col_num_liver(1:num_TimePoint+j)),i);
        end
        p_liver(num_Genotype+j,i) = anova1([temp_list_WT,temp_list_OB],{'WT','OB'},'off');
    end
end
 
%%% q values for those with enough data (WT t0-tX;OB t0-tX;OB/WT tX)(x = 0,2,4,6,8,12,16,24)
%%% ('Method','polynomial','Lambda',0.01:0.01:0.90)or('BHFDR',true)
q_liver = NaN(size(p_liver));
for i = 1:size(p_liver,1)
    q_liver(i,is_EnoughData) = mafdr(p_liver(i,is_EnoughData),'BHFDR',true);
end
 
%%% AUC (log2(data/geometric mean of WT_0 and ob_0))
geomean_time0 = geomean([mean_liver(1,:);mean_liver(num_TimePoint+1,:)],1);
log2normalize_liverWT = log2(mean_liver(1:num_TimePoint,:)./geomean_time0);
AUC_curve_liverWT = log2normalize_liverWT - log2normalize_liverWT(1,:);
log2normalize_liverOB = log2(mean_liver(num_TimePoint+1:end,:)./geomean_time0);
AUC_curve_liverOB = log2normalize_liverOB - log2normalize_liverOB(1,:);
AUC_liverWT = trapz(TimePoint,AUC_curve_liverWT,1);
AUC_liverOB = trapz(TimePoint,AUC_curve_liverOB,1);
 
%%% Responsive & Different (with Enough Data)
is_Responsive(index_WT,:) = (q_liver(index_WT,:)<threshold_q).*(AUC_liverWT>0)*1 + (q_liver(index_WT,:)<threshold_q).*(AUC_liverWT<0)*(-1);
is_Responsive(index_OB,:) = (q_liver(index_OB,:)<threshold_q).*(AUC_liverOB>0)*1 + (q_liver(index_OB,:)<threshold_q).*(AUC_liverOB<0)*(-1);
is_Different(1:num_TimePoint,:) = (q_liver(num_Genotype+1:end,:)<threshold_q).*(FC_liver(num_TimePoint*2+1:num_TimePoint*3,:)>1)*1 + (q_liver(num_Genotype+1:end,:)<threshold_q).*(FC_liver(num_TimePoint*2+1:num_TimePoint*3,:)<1)*(-1);
% Mixed as 2
is_TimeCourseDifferent = false(1,length(is_Different));
is_TimeCourseDifferent = is_TimeCourseDifferent + (sum(logical(is_Different))~=0).*(logical(~sum(is_Different==-1))).*(1);
is_TimeCourseDifferent = is_TimeCourseDifferent + (sum(logical(is_Different))~=0).*(logical(~sum(is_Different==1))).*(-1);
is_TimeCourseDifferent = is_TimeCourseDifferent + (sum(logical(is_Different))~=0).*(logical(sum(is_Different==1)).*logical(sum(is_Different==-1))).*(2);

%%% Calculate T half
% 1=fast,2=medium,3=slow
% {"<8 h","8-16 h",">16 h"}
THalf_liver = NaN(num_Genotype,size(Info_ProteomeName,1));
for i = 1:size(THalf_liver,1)
    for j = 1:size(THalf_liver,2)
        temp_y = mean_liver((i-1)*num_TimePoint+1:num_TimePoint*i,j);
        if is_Responsive(i,j) == 1
            temp_ymax = max(temp_y);
            temp_y0 = mean([temp_y(1),temp_ymax]);
            temp_index_max = find(temp_y == temp_ymax);
            for k = 1:temp_index_max-1
                if (temp_y(k)<=temp_y0)&&(temp_y0<=temp_y(k+1))
                    temp_THalf = GetTHalf(TimePoint(k),temp_y(k),TimePoint(k+1),temp_y(k+1),temp_y0);
                    THalf_liver(i,j) = temp_THalf;
                    break
                end
            end
        elseif is_Responsive(i,j) == -1
            temp_ymin = min(temp_y);
            temp_y0 = mean([temp_y(1),temp_ymin]);
            temp_index_min = find(temp_y == temp_ymin);
            for k = 1:temp_index_min-1
                if (temp_y(k)>=temp_y0)&&(temp_y0>=temp_y(k+1))
                    temp_THalf = GetTHalf(TimePoint(k),temp_y(k),TimePoint(k+1),temp_y(k+1),temp_y0);
                    THalf_liver(i,j) = temp_THalf;
                    break
                end
            end
        end
    end
end
THalfType_liver = (THalf_liver<8)*1+(THalf_liver>=8).*(THalf_liver<=16)*2+(THalf_liver>16)*3;

%% Bimodality Coefficient for THalf
BC_WT = BimodalityCoefficient(THalf_liver(index_WT,logical(is_Responsive(index_WT,:))));
BC_OB = BimodalityCoefficient(THalf_liver(index_OB,logical(is_Responsive(index_OB,:))));
[~,p_AD_WT] = adtest(THalf_liver(index_WT,logical(is_Responsive(index_WT,:))));
[~,p_AD_OB] = adtest(THalf_liver(index_OB,logical(is_Responsive(index_OB,:))));

%% %%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure - Bar Graph (Responsive)
%%% Parameters
bar_color = [1,0.6,0;0.6,0,0.6];
%bar_color = [1,1,1;0.75,0.75,0.75];
name_figure = 'Proteome_Responsive';
%name_figure = 'Proteome_Responsive_B&W';

%%% Draw
bar_data_liver_WT = [sum(is_Responsive(index_WT,:)==1),sum(is_Responsive(index_WT,:)==-1)];
bar_data_liver_ob = [sum(is_Responsive(index_OB,:)==1),sum(is_Responsive(index_OB,:)==-1)];
bar_X = [1,2];
bar_Y = [bar_data_liver_WT;bar_data_liver_ob];
bar_XLabel = ["WT","\it{ob}\rm{/}\it{ob}"];
bar_XLabel_OmicType = "Proteome";
f = figure(1);
hold on;
for i = 1:length(bar_X)
    b = bar(bar_X(i),bar_Y(i,:),'stacked');
    text(bar_X(i),-8,bar_XLabel(i),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',10);
    for j = 1:size(bar_Y(i,:),2)
        temp_Y_pos = sum(bar_Y(i,1:j))-0.5*bar_Y(i,j);
        if bar_Y(i,j) == 0
            continue;
        end
        text(bar_X(i),temp_Y_pos,string(bar_Y(i,j)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',10);
        b(j).FaceColor = bar_color(j,:);
    end
end
for i = 1:num_TissueType
    temp_X_pos = mean([bar_X(i*2),bar_X(i*2-1)]);
    text(temp_X_pos,-35,bar_XLabel_OmicType(i),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',15)
end
temp_loc = 0.5*max([sum(bar_data_liver_WT),sum(bar_data_liver_ob)]);
text(-1.2,temp_loc,'Counts','HorizontalAlignment','center','VerticalAlignment','middle','Rotation',90,'FontSize',15);
set(gca,'XColor','none');
legend(b,{'Increase','Decrease'});
%ylim([0,100]);
f.Position = [50 50 400 400];

%% Save
exportgraphics(f,['./Fig/',name_figure,'.pdf'],'ContentType','vector');

%% Figure - Bar Graph (T Half)
%%% Parameters
threshold_THalf_1 = 8;
threshold_THalf_2 = 16;
%bar_color = [0.8,1,0.8;1,1,0.8;1,0.8,0.8];
bar_color = [1,1,1;0.75,0.75,0.75;0.5,0.5,0.5];
%name_figure = 'Proteome_THalf';
name_figure = 'Proteome_THalf_B&W';

%%% Draw
bar_data_liver_WT = [sum(THalf_liver(1,:)<threshold_THalf_1),sum((THalf_liver(1,:)>=threshold_THalf_1).*(THalf_liver(1,:)<=threshold_THalf_2)),sum(THalf_liver(1,:)>threshold_THalf_2)];
bar_data_liver_ob = [sum(THalf_liver(2,:)<threshold_THalf_1),sum((THalf_liver(2,:)>=threshold_THalf_1).*(THalf_liver(2,:)<=threshold_THalf_2)),sum(THalf_liver(2,:)>threshold_THalf_2)];
bar_X = [1,2];
bar_Y = [bar_data_liver_WT;bar_data_liver_ob];
bar_XLabel = ["WT","\it{ob}\rm{/}\it{ob}"];
bar_XLabel_OmicType = "Proteome";
f = figure(1);
hold on;
for i = 1:length(bar_X)
    b = bar(bar_X(i),bar_Y(i,:),'stacked');
    text(bar_X(i),-8,bar_XLabel(i),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',10);
    for j = 1:size(bar_Y(i,:),2)
        temp_Y_pos = sum(bar_Y(i,1:j))-0.5*bar_Y(i,j);
        if bar_Y(i,j) == 0
            continue;
        end
        text(bar_X(i),temp_Y_pos,string(bar_Y(i,j)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',10);
        b(j).FaceColor = bar_color(j,:);
    end
end
for i = 1:num_TissueType
    temp_X_pos = mean([bar_X(i*2),bar_X(i*2-1)]);
    text(temp_X_pos,-35,bar_XLabel_OmicType(i),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',15)
end
temp_loc = 0.5*max([sum(bar_data_liver_WT),sum(bar_data_liver_ob)]);
text(-1.2,temp_loc,'Counts','HorizontalAlignment','center','VerticalAlignment','middle','Rotation',90,'FontSize',15);
set(gca,'XColor','none');
legend(b,{['< ',num2str(threshold_THalf_1),' h'],[num2str(threshold_THalf_1),'-',num2str(threshold_THalf_2),' h'],['> ',num2str(threshold_THalf_2),' h']});
%ylim([0,100]);
f.Position = [50 50 400 400];

%% Save
exportgraphics(f,['./Fig/',name_figure,'.pdf'],'ContentType','vector');

%% %%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure - Heat map
%% Parameters
color_low = [0.4940,0.1840,0.5560];
color_zero = [1,1,1];
color_high = [0.9290,0.6940,0.1250];
color_WT = [0,0,1];
color_OB = [1,0,0];
color_fast = [0.8,1,0.8];
color_medium = [1,1,0.8];
color_slow = [1,0.8,0.8];

color_mix = [1,0.8,1];
color_OB_high = [1,0.8,0.8];
color_WT_high = [0.8,0.8,1];

prefix_common = '\color[rgb]{0,1,0}';
prefix_WT = '\color[rgb]{0,0,1}';
prefix_OB = '\color[rgb]{1,0,0}';
prefix_opposite = '\color[rgb]{1,0,1}';
prefix_unresp = '\color[rgb]{0,0,0}';
color_Scale = zeros(200,length(color_zero));
for i = 1:length(color_zero)
    temp_list = [linspace(color_low(i),color_zero(i)),linspace(color_zero(i),color_high(i))];
    color_Scale(:,i) = temp_list';
end
%cutoff = 0.5;
temp_num_Cluster = 3;
name_figure_1 = 'Proteome_liver';
name_figure_2 = 'Proteome_MeanTimeCourseInOne_liver';
%name_figure_suffix = ['_',num2str(cutoff),'MaxLinkage'];
name_figure_suffix = ['_',num2str(temp_num_Cluster),'_Clusters'];

%% Proteome_liver_WT&OB
%%% Dendrogram
f1 = figure(1);
subplot(1,4,1);
temp_include = logical(is_EnoughData.*(sum(logical(is_Responsive),1)>0));
temp_data = [log2normalize_liverWT(:,temp_include);log2normalize_liverOB(:,temp_include)];
[temp_data_sorted,temp_index_row,temp_tree] = SortRowByCluster(temp_data','ward','euclidean');
temp_index_row_dendrogram = flip(temp_index_row);
temp_YLabel = cell(1,length(temp_data));
for i = 1:length(temp_data)
    temp_YLabel{i} = '';
end

cutoff = median([temp_tree(end-temp_num_Cluster+1,3),temp_tree(end-temp_num_Cluster+2,3)]);
[H,~,~] = dendrogram(temp_tree,0,'Orientation','left','Reorder',temp_index_row_dendrogram,'ColorThreshold',cutoff,'Labels',temp_YLabel); 

%[H,~,~] = dendrogram(temp_tree,0,'Orientation','left','Reorder',temp_index_row_dendrogram,'ColorThreshold',cutoff*max(temp_tree(:,3)),'Labels',temp_YLabel); 
set(gca, 'XColor','none', 'YColor','none','Color','none');
%%% Node - Cluster - Color [*,*,*]
temp_Info_Cluster = zeros(length(temp_index_row),5);
temp_Info_Cluster(:,1) = temp_index_row';
for i = 1:size(temp_Info_Cluster,1)
    temp_temp_node = temp_Info_Cluster(i,1);
    [temp_temp_row,~] = find(temp_tree(:,1:2) == temp_temp_node);
    temp_Info_Cluster(i,3:5) = H(temp_temp_row).Color;
end
for i = 1:size(temp_Info_Cluster,1)
    if i == 1
        temp_color_compare = temp_Info_Cluster(1,3:5);
        temp_Info_Cluster(1,2) = 1;
        continue
    elseif prod(temp_color_compare == temp_Info_Cluster(i,3:5)) == 1
        temp_Info_Cluster(i,2) = temp_Info_Cluster(i-1,2);
    else
        temp_Info_Cluster(i,2) = temp_Info_Cluster(i-1,2)+1;
    end
    temp_color_compare = temp_Info_Cluster(i,3:5);
end
%%% Text
text(0.1,0.5,[num2str(sum(temp_include)),' proteins'],'FontSize',10);
text(0.1,0.3,"WT",'FontSize',10);
text(0.5,0.3,"\it{ob}\rm{/}\it{ob}",'FontSize',10);
text(0.3,0.2,"(h)",'FontSize',5);
%%% Change Dendrogram Color
temp_Colormap = [0,0,0;parula(8)];
%temp_Colormap = temp_Colormap([5,6,8,4,7,3],:);
temp_Colormap = temp_Colormap([1,4,6,7],:);
ChangeDendrogramColor(H,temp_Colormap);

%%% Heat map 1
subplot(1,4,2);
temp_XLabel = {'0','2','4','6','8','12','16','24','0','2','4','6','8','12','16','24'};
temp_YLabel = Info_ProteomeName(temp_include,1);
for i = 1:length(temp_YLabel)
    temp_YLabel{i} = '';
end
h1 = heatmap(temp_data_sorted,'ColorLimits',[-2,2],'Colormap',color_Scale,'CellLabelColor','none');
h1.XDisplayLabels = temp_XLabel;
h1.YDisplayLabels = temp_YLabel';
h1.GridVisible = 'off';
h1.FontSize = 5;

%%% Heat map 2
subplot(1,4,3);
temp_data = is_Responsive(:,temp_include);
temp_data = temp_data';
temp_data_sorted = temp_data(temp_index_row,:);
h2 = heatmap(temp_data_sorted,'ColorLimits',[-1,1],'Colormap',[color_low;color_zero;color_high],'CellLabelColor','none');
h2.ColorbarVisible = 'off';
h2.GridVisible = 'off';
h2.XDisplayLabels = ["WT","\it{ob}\rm{/}\it{ob}"];
h2.YDisplayLabels = temp_YLabel;
h2.FontSize = 5;

%% Heat map 3
%{
subplot(1,4,4);
temp_data = is_TimeCourseDifferent(:,temp_include);
temp_data = temp_data';
temp_data_sorted = temp_data(temp_index_row,:);
h3 = heatmap(temp_data_sorted,'ColorLimits',[-1,2],'Colormap',[color_WT_high;color_zero;color_OB_high;color_mix],'CellLabelColor','none');
h3.ColorbarVisible = 'off';
h3.GridVisible = 'off';
h3.XDisplayLabels = "Dominance";
for i = 1:length(temp_YLabel)
    temp_YLabel{i} = '';
end
h3.YDisplayLabels = temp_YLabel;
h3.FontSize = 5;

f1.Position = [10 10 1000 700];
%}

%% Heat map 3
subplot(1,4,4);
temp_data = THalfType_liver(:,temp_include);
temp_data = temp_data';
temp_data_sorted = temp_data(temp_index_row,:);
h3 = heatmap(temp_data_sorted,'ColorLimits',[0,3],'Colormap',[color_zero;color_fast;color_medium;color_slow],'CellLabelColor','none');
h3.ColorbarVisible = 'off';
h3.GridVisible = 'off';
h3.XDisplayLabels = ["WT","\it{ob}\rm{/}\it{ob}"];
for i = 1:length(temp_YLabel)
    temp_YLabel{i} = '';
end
h3.YDisplayLabels = temp_YLabel;
h3.FontSize = 5;

%f1.Position = [10 10 1000 700];

%% Mean Time Course - All Tissue in One
f2 = figure(2);
temp_x = TimePoint;
temp_include = logical(is_EnoughData.*(sum(logical(is_Responsive),1)>0));
temp_data = [log2normalize_liverWT(:,temp_include);log2normalize_liverOB(:,temp_include)];
temp_data = temp_data';
temp_Cluster_list = unique(temp_Info_Cluster(:,2),'stable');
for i = 1:temp_Cluster_list(end)
    temp_temp_node_list = temp_Info_Cluster(temp_Info_Cluster(:,2)==temp_Cluster_list(i),1);
    temp_data_cluster = temp_data(temp_temp_node_list,:);
    temp_data_cluster_mean = mean(temp_data_cluster,1);
    temp_data_cluster_error = std(temp_data_cluster,0,1);
    subplot(temp_Cluster_list(end),1,i);
    xticks(temp_x);
    xlabel(["Time (h)";strcat('Cluster ',num2str(i))]);
    ylabel("log2FC");
    xlim([0,24]);
    hold on;
    % WT
    errorbar(temp_x,temp_data_cluster_mean(1:num_TimePoint),temp_data_cluster_error(1:num_TimePoint),'Color',color_WT,'LineWidth',1);
    % OB
    errorbar(temp_x,temp_data_cluster_mean(num_TimePoint+1:num_TimePoint*2),temp_data_cluster_error(num_TimePoint+1:num_TimePoint*2),'Color',color_OB,'LineWidth',1);
end
legend(["WT","\it{ob}\rm{/}\it{ob}"]);
f2.Position = [0 0 200 700];

%% Save
exportgraphics(f1,['./Fig/',name_figure_1,name_figure_suffix,'.pdf'],'ContentType','vector');
%%
exportgraphics(f2,['./Fig/',name_figure_2,name_figure_suffix,'.pdf'],'ContentType','vector');

%% %%%%%%%%%%%%%%%%%%%%%%%%%
%% Function
function logic_EnoughData = CheckEnoughData(col_var, col_num)
logic_EnoughData = 1;
temp_sum = 0;
for i = 1:length(col_num)
    temp_data = col_var(temp_sum+1:sum(col_num(1:i)));
    temp_threshold = col_num(i)/2;
    logic_EnoughData = logic_EnoughData*(sum(isnan(temp_data)) < temp_threshold);
    temp_sum = temp_sum + col_num(i);
end
logic_EnoughData = logical(logic_EnoughData);
end

function result_mean = GetMean(col_var, col_num)
result_mean = NaN(length(col_num),1);
temp_sum = 0;
for i = 1:length(col_num)
    temp_data = col_var(temp_sum+1:sum(col_num(1:i)));
    result_mean(i) = mean(temp_data,'omitnan');
    temp_sum = temp_sum + col_num(i);
end
end

function [data_sorted,index_row,tree] = SortRowByCluster(data,method,metric)
if nargin < 2
   method = 'average';
   metric = 'correlation';
end
tree = linkage(data,method,metric);
D = pdist(data,metric);
index_row = optimalleaforder(tree,D);
data_sorted = data(index_row,:);
end

function THalf = GetTHalf(x1,y1,x2,y2,y0)
if y0 == y1
    THalf = x1;
elseif y0 == y2
    THalf = x2;
else
    THalf = x2 + (x2-x1)*(y0-y2)/(y2-y1);
end
end

function ChangeDendrogramColor(H,Colormap)
DendrogramColor = cell2mat(get(H,'Color'));
DendrogramColor_unique = unique(DendrogramColor,'row');
for i = 1:size(H,1)
    H(i).Color = Colormap(logical(prod(DendrogramColor_unique == H(i).Color,2)),:);
end
end

function BC = BimodalityCoefficient(data)
g = skewness(data,0);
k = kurtosis(data,0);
n = length(data);
BC = (g.^2+1)/(k+3*(n-1).^2/(n-2)/(n-3));
end
