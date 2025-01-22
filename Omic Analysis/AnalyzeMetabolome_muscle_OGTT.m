%% Load files
load('./Resource/MetabolomeData_OGTT_muscle&liver.mat');

%% Parameters
num_Genotype = 2;
index_WT = 1;
index_OB = 2;
num_TissueType = 1;
index_muscle = 1;
num_TimePoint = 5;
TimePoint = [0,20,60,120,240];
col_num_muscle = [5,5,5,5,5,5,5,5,5,5];
num_Replicate = 5;
Group = {'WT0','WT0','WT0','WT0','WT0','WT20','WT20','WT20','WT20','WT20','WT60','WT60','WT60','WT60','WT60','WT120','WT120','WT120','WT120','WT120','WT240','WT240','WT240','WT240','WT240','OB0','OB0','OB0','OB0','OB0','OB20','OB20','OB20','OB20','OB20','OB60','OB60','OB60','OB60','OB60','OB120','OB120','OB120','OB120','OB120','OB240','OB240','OB240','OB240','OB240'};
threshold_q = 0.1;

%% Check Enough Data (Exist in >= half replicates)
is_EnoughData = false(num_TissueType,size(Info_MetabolomeName,1));
for i = 1:size(is_EnoughData,2)
    logic_EnoughData = CheckEnoughData(Metabolome_muscle(:,i), col_num_muscle);
    is_EnoughData(index_muscle,i) = logic_EnoughData;
end

%% Calculation
%%% Mean
mean_muscle = NaN(length(col_num_muscle),size(Info_MetabolomeName,1));
for i = 1:size(mean_muscle,2)
    mean_muscle(:,i) = GetMean(Metabolome_muscle(:,i),col_num_muscle);
end

%%% FC & log2FC (WT tX/t0, OB tX/t0, OB/WT tX)(x = 0,2,4,6,8,12,16,24)
FC_muscle = [mean_muscle(1:num_TimePoint,:)./mean_muscle(1,:);mean_muscle(num_TimePoint+1:num_TimePoint*2,:)./mean_muscle(num_TimePoint+1,:);mean_muscle(num_TimePoint+1:num_TimePoint*2,:)./mean_muscle(1:num_TimePoint,:)];
log2FC_muscle = log2(FC_muscle);

%%% p values (One-way ANOVA)(WT;OB;OB/WT tX)
p_muscle = NaN(num_Genotype,size(Info_MetabolomeName,1));
temp_num = sum(col_num_muscle(1:num_TimePoint));
% muscle - WT
for i = 1:size(Info_MetabolomeName,1)
    p_muscle(index_WT,i) = anova1(Metabolome_muscle(1:temp_num,i),Group(1:temp_num),'off');
end
% muscle - OB
for i = 1:size(Info_MetabolomeName,1)
    p_muscle(index_OB,i) = anova1(Metabolome_muscle(temp_num+1:end,i),Group(temp_num+1:end),'off');
end
% muscle - OB/WT
for i = 1:size(Info_MetabolomeName,1)
    for j = 1:num_TimePoint
        if j == 1
            temp_list_WT = Metabolome_muscle(1:col_num_muscle(j),i);
            temp_list_OB = Metabolome_muscle(temp_num+1:temp_num+col_num_muscle(num_TimePoint+j),i);
        else
            temp_list_WT = Metabolome_muscle(sum(col_num_muscle(1:j-1))+1:sum(col_num_muscle(1:j)),i);
            temp_list_OB = Metabolome_muscle(sum(col_num_muscle(1:num_TimePoint+j-1))+1:sum(col_num_muscle(1:num_TimePoint+j)),i);
        end
        p_muscle(num_Genotype+j,i) = anova1([temp_list_WT,temp_list_OB],{'WT','OB'},'off');
    end
end

%%% q values for those with enough data (WT t0-tX;OB t0-tX;OB/WT tX)(x = 0,2,4,6,8,12,16,24)
%%% ('Method','polynomial','Lambda',0.01:0.01:0.90)or('BHFDR',true)
q_muscle = NaN(size(p_muscle));
for i = 1:size(p_muscle,1)
    q_muscle(i,is_EnoughData) = mafdr(p_muscle(i,is_EnoughData),'BHFDR',true);
end

%%% AUC (log2(data/geometric mean of WT_0 and ob_0))
geomean_time0 = geomean([mean_muscle(1,:);mean_muscle(num_TimePoint+1,:)],1);
log2normalize_muscleWT = log2(mean_muscle(1:num_TimePoint,:)./geomean_time0);
AUC_curve_muscleWT = log2normalize_muscleWT - log2normalize_muscleWT(1,:);
log2normalize_muscleOB = log2(mean_muscle(num_TimePoint+1:end,:)./geomean_time0);
AUC_curve_muscleOB = log2normalize_muscleOB - log2normalize_muscleOB(1,:);
AUC_muscleWT = trapz(TimePoint,AUC_curve_muscleWT,1);
AUC_muscleOB = trapz(TimePoint,AUC_curve_muscleOB,1);

%%% Responsive & Different (with Enough Data)
is_Responsive(index_WT,:) = (q_muscle(index_WT,:)<threshold_q).*(AUC_muscleWT>0)*1 + (q_muscle(index_WT,:)<threshold_q).*(AUC_muscleWT<0)*(-1);
is_Responsive(index_OB,:) = (q_muscle(index_OB,:)<threshold_q).*(AUC_muscleOB>0)*1 + (q_muscle(index_OB,:)<threshold_q).*(AUC_muscleOB<0)*(-1);
is_Different(1:num_TimePoint,:) = (q_muscle(num_Genotype+1:end,:)<threshold_q).*(FC_muscle(num_TimePoint*2+1:num_TimePoint*3,:)>1)*1 + (q_muscle(num_Genotype+1:end,:)<threshold_q).*(FC_muscle(num_TimePoint*2+1:num_TimePoint*3,:)<1)*(-1);
% Mixed as 2
is_TimeCourseDifferent = false(1,size(is_Different,2));
is_TimeCourseDifferent = is_TimeCourseDifferent + (sum(logical(is_Different))~=0).*(logical(~sum(is_Different==-1))).*(1);
is_TimeCourseDifferent = is_TimeCourseDifferent + (sum(logical(is_Different))~=0).*(logical(~sum(is_Different==1))).*(-1);
is_TimeCourseDifferent = is_TimeCourseDifferent + (sum(logical(is_Different))~=0).*(logical(sum(is_Different==1)).*logical(sum(is_Different==-1))).*(2);

%%% Calculate T half
% 1=fast,2=medium,3=slow
% {"<20 min","20-60 min",">60 min"}
THalf_muscle = NaN(num_Genotype,size(Info_MetabolomeName,1));
for i = 1:size(THalf_muscle,1)
    for j = 1:size(THalf_muscle,2)
        temp_y = mean_muscle((i-1)*num_TimePoint+1:num_TimePoint*i,j);
        if is_Responsive(i,j) == 1
            temp_ymax = max(temp_y);
            temp_y0 = mean([temp_y(1),temp_ymax]);
            temp_index_max = find(temp_y == temp_ymax);
            for k = 1:temp_index_max-1
                if (temp_y(k)<=temp_y0)&&(temp_y0<=temp_y(k+1))
                    temp_THalf = GetTHalf(TimePoint(k),temp_y(k),TimePoint(k+1),temp_y(k+1),temp_y0);
                    THalf_muscle(i,j) = temp_THalf;
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
                    THalf_muscle(i,j) = temp_THalf;
                    break
                end
            end
        end
    end
end
THalfType_muscle = (THalf_muscle<20)*1+(THalf_muscle>=20).*(THalf_muscle<=60)*2+(THalf_muscle>60)*3;

%% Energy Charge & AMP/ATP & ADP/ATP
index_ATP = 1;
index_ADP = 6;
index_AMP = 11;
% WT
ATP = Metabolome_muscle(1:num_Replicate*num_TimePoint,index_ATP);
ADP = Metabolome_muscle(1:num_Replicate*num_TimePoint,index_ADP);
AMP = Metabolome_muscle(1:num_Replicate*num_TimePoint,index_AMP);
EC_WT = (ATP + 0.5.*ADP)./(ATP + ADP + AMP);
AMP2ATP_WT = AMP./ATP;
ADP2ATP_WT = ADP./ATP;
p_EC_WT = anova1(EC_WT,Group(1:num_Replicate*num_TimePoint),'off');
p_AMP2ATP_WT = anova1(AMP2ATP_WT,Group(1:num_Replicate*num_TimePoint),'off');
p_ADP2ATP_WT = anova1(ADP2ATP_WT,Group(1:num_Replicate*num_TimePoint),'off');
% OB
ATP = Metabolome_muscle(num_Replicate*num_TimePoint+1:end,index_ATP);
ADP = Metabolome_muscle(num_Replicate*num_TimePoint+1:end,index_ADP);
AMP = Metabolome_muscle(num_Replicate*num_TimePoint+1:end,index_AMP);
EC_OB = (ATP + 0.5.*ADP)./(ATP + ADP + AMP);
AMP2ATP_OB = AMP./ATP;
ADP2ATP_OB = ADP./ATP;
p_EC_OB = anova1(EC_OB,Group(num_Replicate*num_TimePoint+1:end),'off');
p_AMP2ATP_OB = anova1(AMP2ATP_OB,Group(num_Replicate*num_TimePoint+1:end),'off');
p_ADP2ATP_OB = anova1(ADP2ATP_OB,Group(num_Replicate*num_TimePoint+1:end),'off');
% OB/WT
p_EC_Diff = NaN(num_TimePoint,1);
for j = 1:num_TimePoint
    if j == 1
        temp_list_WT = EC_WT(1:col_num_muscle(j));
        temp_list_OB = EC_OB(1:col_num_muscle(j));
    else
        temp_list_WT = EC_WT(sum(col_num_muscle(1:j-1))+1:sum(col_num_muscle(1:j)));
        temp_list_OB = EC_OB(sum(col_num_muscle(1:j-1))+1:sum(col_num_muscle(1:j)));
    end
    p_EC_Diff(j) = anova1([temp_list_WT,temp_list_OB],{'WT','OB'},'off');
end
p_AMP2ATP_Diff = NaN(num_TimePoint,1);
for j = 1:num_TimePoint
    if j == 1
        temp_list_WT = AMP2ATP_WT(1:col_num_muscle(j));
        temp_list_OB = AMP2ATP_OB(1:col_num_muscle(j));
    else
        temp_list_WT = AMP2ATP_WT(sum(col_num_muscle(1:j-1))+1:sum(col_num_muscle(1:j)));
        temp_list_OB = AMP2ATP_OB(sum(col_num_muscle(1:j-1))+1:sum(col_num_muscle(1:j)));
    end
    p_AMP2ATP_Diff(j) = anova1([temp_list_WT,temp_list_OB],{'WT','OB'},'off');
end
p_ADP2ATP_Diff = NaN(num_TimePoint,1);
for j = 1:num_TimePoint
    if j == 1
        temp_list_WT = ADP2ATP_WT(1:col_num_muscle(j));
        temp_list_OB = ADP2ATP_OB(1:col_num_muscle(j));
    else
        temp_list_WT = ADP2ATP_WT(sum(col_num_muscle(1:j-1))+1:sum(col_num_muscle(1:j)));
        temp_list_OB = ADP2ATP_OB(sum(col_num_muscle(1:j-1))+1:sum(col_num_muscle(1:j)));
    end
    p_ADP2ATP_Diff(j) = anova1([temp_list_WT,temp_list_OB],{'WT','OB'},'off');
end

%% Figure - Energy Charge & AMP/ATP & ADP/ATP
%{
%%% EC
EC_WT_matrix = reshape(EC_WT,num_Replicate,num_TimePoint);
EC_OB_matrix = reshape(EC_OB,num_Replicate,num_TimePoint);
EC_WT_mean = mean(EC_WT_matrix,1);
EC_OB_mean = mean(EC_OB_matrix,1);
EC_WT_error = std(EC_WT_matrix,0,1);
EC_OB_error = std(EC_OB_matrix,0,1);

f1 = figure(1);
hold on;
for i = 1:size(EC_WT_matrix,1)
    scatter(TimePoint,EC_WT_matrix(i,:),[],[0.8 0.8 1],'filled');
end
for i = 1:size(EC_OB_matrix,1)
    scatter(TimePoint,EC_OB_matrix(i,:),[],[1 0.8 0.8],'filled');
end
p1 = errorbar(TimePoint,EC_WT_mean,EC_WT_error,'Color',[0 0 1],'LineWidth',1);
p2 = errorbar(TimePoint,EC_OB_mean,EC_OB_error,'Color',[1 0 0],'LineWidth',1);
text(8,0.93,'*','FontSize',20,'Color','k','HorizontalAlignment','center');          % Manual

xlim([0,24]);
xticks(TimePoint);
xlabel("Time Course (h)");
ylim([0.8,1]);                                                                      % Manual
ylabel("Energy Charge");
legend([p1,p2],["WT","\it{ob}\rm{/}\it{ob}"]);
f1.Position = [100 500 220 200];
%}
%%% AMP/ATP
AMP2ATP_WT_matrix = reshape(AMP2ATP_WT,num_Replicate,num_TimePoint);
AMP2ATP_OB_matrix = reshape(AMP2ATP_OB,num_Replicate,num_TimePoint);
AMP2ATP_WT_mean = mean(AMP2ATP_WT_matrix,1);
AMP2ATP_OB_mean = mean(AMP2ATP_OB_matrix,1);
AMP2ATP_WT_error = std(AMP2ATP_WT_matrix,0,1);
AMP2ATP_OB_error = std(AMP2ATP_OB_matrix,0,1);

f2 = figure(2);
hold on;
for i = 1:size(AMP2ATP_WT_matrix,1)
    scatter(TimePoint,AMP2ATP_WT_matrix(i,:),[],[0.8 0.8 1],'filled');
end
for i = 1:size(AMP2ATP_OB_matrix,1)
    scatter(TimePoint,AMP2ATP_OB_matrix(i,:),[],[1 0.8 0.8],'filled');
end
p1 = errorbar(TimePoint,AMP2ATP_WT_mean,AMP2ATP_WT_error,'Color',[0 0 1],'LineWidth',1);
p2 = errorbar(TimePoint,AMP2ATP_OB_mean,AMP2ATP_OB_error,'Color',[1 0 0],'LineWidth',1);
%text(8,0.07,'*','FontSize',20,'Color','k','HorizontalAlignment','center');                      % Manual
%text(5,0.09,'\bf{WT}','FontSize',15,'Color','b','HorizontalAlignment','center');             % Manual
%text(5.4,0.095,'\bf{WT}','FontSize',12,'Color','b','HorizontalAlignment','center');             % Manual
%text(5.4,0.085,'\bf{increase}','FontSize',12,'Color','b','HorizontalAlignment','center');       % Manual

xlim([0,240]);
xticks(TimePoint);
xlabel("Time Course (min)");
ylim([0,0.04]);                                                                                  % Manual
ylabel("AMP/ATP");
legend([p1,p2],["WT","\it{ob}\rm{/}\it{ob}"]);
f2.Position = [500 500 220 200];

%{
%%% ADP/ATP
ADP2ATP_WT_matrix = reshape(ADP2ATP_WT,num_Replicate,num_TimePoint);
ADP2ATP_OB_matrix = reshape(ADP2ATP_OB,num_Replicate,num_TimePoint);
ADP2ATP_WT_mean = mean(ADP2ATP_WT_matrix,1);
ADP2ATP_OB_mean = mean(ADP2ATP_OB_matrix,1);
ADP2ATP_WT_error = std(ADP2ATP_WT_matrix,0,1);
ADP2ATP_OB_error = std(ADP2ATP_OB_matrix,0,1);

f3 = figure(3);
hold on;
for i = 1:size(ADP2ATP_WT_matrix,1)
    scatter(TimePoint,ADP2ATP_WT_matrix(i,:),[],[0.8 0.8 1],'filled');
end
for i = 1:size(ADP2ATP_OB_matrix,1)
    scatter(TimePoint,ADP2ATP_OB_matrix(i,:),[],[1 0.8 0.8],'filled');
end
p1 = errorbar(TimePoint,ADP2ATP_WT_mean,ADP2ATP_WT_error,'Color',[0 0 1],'LineWidth',1);
p2 = errorbar(TimePoint,ADP2ATP_OB_mean,ADP2ATP_OB_error,'Color',[1 0 0],'LineWidth',1);
text(8,0.3,'*','FontSize',20,'Color','k','HorizontalAlignment','center');          % Manual

xlim([0,24]);
xticks(TimePoint);
xlabel("Time Course (h)");
ylim([0.1,0.4]);                                                                      % Manual
ylabel("ADP/ATP");
legend([p1,p2],["WT","\it{ob}\rm{/}\it{ob}"]);
f3.Position = [900 500 220 200];
%}

%% Save
exportgraphics(f2,['./Fig/','OGTT_muscle_AMP2ATP','.pdf'],'ContentType','vector');

%% %%%%%%%%%%%%%%%%%%%%%%%%%
%{
%% Figure - Bar Graph (Responsive)
%%% Parameters
bar_color = [1,0.6,0;0.6,0,0.6];
%bar_color = [1,1,1;0.75,0.75,0.75];
name_figure = 'Metabolome_Responsive';
%name_figure = 'Metabolome_Responsive_B&W';

%%% Draw
bar_data_muscle_WT = [sum(is_Responsive(index_WT,:)==1),sum(is_Responsive(index_WT,:)==-1)];
bar_data_muscle_ob = [sum(is_Responsive(index_OB,:)==1),sum(is_Responsive(index_OB,:)==-1)];
bar_X = [1,2];
bar_Y = [bar_data_muscle_WT;bar_data_muscle_ob];
bar_XLabel = ["WT","\it{ob}\rm{/}\it{ob}"];
bar_XLabel_TissueType = "Muscle";
f = figure(1);
hold on;
for i = 1:length(bar_X)
    b = bar(bar_X(i),bar_Y(i,:),'stacked');
    text(bar_X(i),-2,bar_XLabel(i),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',10);
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
    text(temp_X_pos,-7,bar_XLabel_TissueType,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',15)
end
temp_loc = 0.5*max([sum(bar_data_muscle_WT),sum(bar_data_muscle_ob)]);
text(-1.2,temp_loc,'Counts','HorizontalAlignment','center','VerticalAlignment','middle','Rotation',90,'FontSize',15);
set(gca,'XColor','none');
legend(b,{'Increase','Decrease'});
ylim([0,100]);
f.Position = [50 50 400 400];

%% Save
exportgraphics(f,['./Fig/',name_figure,'.pdf'],'ContentType','vector');

%% Figure - Bar Graph (T Half)
%%% Parameters
threshold_THalf_1 = 8;
threshold_THalf_2 = 16;
bar_color = [0.8,1,0.8;1,1,0.8;1,0.8,0.8];
%bar_color = [1,1,1;0.75,0.75,0.75;0.5,0.5,0.5];
name_figure = 'Metabolome_THalf';
%name_figure = 'Metabolome_THalf_B&W';

%%% Draw
bar_data_muscle_WT = [sum(THalf_muscle(index_WT,:)<threshold_THalf_1),sum((THalf_muscle(index_WT,:)>=threshold_THalf_1).*(THalf_muscle(index_WT,:)<=threshold_THalf_2)),sum(THalf_muscle(index_WT,:)>threshold_THalf_2)];
bar_data_muscle_ob = [sum(THalf_muscle(index_OB,:)<threshold_THalf_1),sum((THalf_muscle(index_OB,:)>=threshold_THalf_1).*(THalf_muscle(index_OB,:)<=threshold_THalf_2)),sum(THalf_muscle(index_OB,:)>threshold_THalf_2)];
bar_X = [1,2];
bar_Y = [bar_data_muscle_WT;bar_data_muscle_ob];
bar_XLabel = ["WT","\it{ob}\rm{/}\it{ob}"];
bar_XLabel_TissueType = "Muscle";
f = figure(1);
hold on;
for i = 1:length(bar_X)
    b = bar(bar_X(i),bar_Y(i,:),'stacked');
    text(bar_X(i),-2,bar_XLabel(i),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',10);
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
    text(temp_X_pos,-7,bar_XLabel_TissueType,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',15)
end
temp_loc = 0.5*max([sum(bar_data_muscle_WT),sum(bar_data_muscle_ob)]);
text(-1.2,temp_loc,'Counts','HorizontalAlignment','center','VerticalAlignment','middle','Rotation',90,'FontSize',15);
set(gca,'XColor','none');
legend(b,{['< ',num2str(threshold_THalf_1),' h'],[num2str(threshold_THalf_1),'-',num2str(threshold_THalf_2),' h'],['> ',num2str(threshold_THalf_2),' h']});
ylim([0,100]);
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
cutoff = 0.5;
name_figure_1 = 'Metabolome_muscle';
name_figure_2 = 'Metabolome_MeanTimeCourseInOne_Muscle';
name_figure_suffix = ['_',num2str(cutoff),'MaxLinkage'];

%% Metabolome_muscle_WT&OB
%%% Dendrogram
f1 = figure(1);
subplot(1,4,1);
temp_include = is_EnoughData;
temp_data = [log2normalize_muscleWT(:,temp_include);log2normalize_muscleOB(:,temp_include)];
[temp_data_sorted,temp_index_row,temp_tree] = SortRowByCluster(temp_data','ward','euclidean');
temp_index_row_dendrogram = flip(temp_index_row);
temp_YLabel = cell(1,length(temp_data));
for i = 1:length(temp_data)
    temp_YLabel{i} = '';
end
[H,~,~] = dendrogram(temp_tree,0,'Orientation','left','Reorder',temp_index_row_dendrogram,'ColorThreshold',cutoff*max(temp_tree(:,3)),'Labels',temp_YLabel); 
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
text(-5,sum(temp_include)+2,[num2str(sum(temp_include)),' metabolites'],'FontSize',10);
text(-5,-3,"WT",'FontSize',10);
text(-10,-3,"\it{ob}\rm{/}\it{ob}",'FontSize',10);
text(-25,-2,"(h)",'FontSize',5);
%%% Change Dendrogram Color
temp_Colormap = copper(8);
%temp_Colormap = temp_Colormap([2,8,6,4,5],:);
temp_Colormap = temp_Colormap([1,5,8,6,4],:);

ChangeDendrogramColor(H,temp_Colormap);

%%% Heat map 1
subplot(1,4,2);
temp_XLabel = {'0','2','4','6','8','12','16','24','0','2','4','6','8','12','16','24'};
temp_YLabel = Info_MetabolomeName(:,2);
for i = 1:length(temp_YLabel)
    temp_list = is_Responsive(1:num_Genotype,i);
    if (temp_list(1)~=0)&&(temp_list(1)==temp_list(2))
        temp_prefix = prefix_common;
    elseif (temp_list(1)~=0)&&(temp_list(2)==0)
        temp_prefix = prefix_WT;
    elseif (temp_list(1)==0)&&(temp_list(2)~=0)
        temp_prefix = prefix_OB;
    elseif (temp_list(1)~=0)&&(temp_list(2)~=0)&&(temp_list(1)~=temp_list(2))
        temp_prefix = prefix_opposite;
    else
        temp_prefix = prefix_unresp;
    end
    temp_YLabel{i} = [temp_prefix,temp_YLabel{i}];
end
temp_YLabel = temp_YLabel(temp_include);
temp_YLabel = temp_YLabel(temp_index_row);

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
for i = 1:length(temp_YLabel)
    temp_YLabel{i} = '';
end
h2.YDisplayLabels = temp_YLabel;
h2.FontSize = 5;

%%% Heat map 3
subplot(1,4,4);
temp_data = THalfType_muscle(:,temp_include);
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

f1.Position = [10 10 1000 900];

%% Mean Time Course - All Genotype in One
f2 = figure(2);
temp_x = TimePoint;
temp_include = is_EnoughData;
temp_data = [log2normalize_muscleWT(:,temp_include);log2normalize_muscleOB(:,temp_include)];
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
%}
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
%{
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
%}
function THalf = GetTHalf(x1,y1,x2,y2,y0)
if y0 == y1
    THalf = x1;
elseif y0 == y2
    THalf = x2;
else
    THalf = x2 + (x2-x1)*(y0-y2)/(y2-y1);
end
end
%{
function ChangeDendrogramColor(H,Colormap)
DendrogramColor = cell2mat(get(H,'Color'));
DendrogramColor_unique = unique(DendrogramColor,'row');
for i = 1:size(H,1)
    H(i).Color = Colormap(logical(prod(DendrogramColor_unique == H(i).Color,2)),:);
end
end
%}