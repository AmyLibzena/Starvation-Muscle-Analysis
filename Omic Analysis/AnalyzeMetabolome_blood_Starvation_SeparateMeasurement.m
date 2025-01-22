%% Load files
load('./Resource/BloodInsulinGlucose.mat');
% Insulin: WT = pg/ml; OB = ng/ml
% Glucose: mg/dl

%% Convert OB insulin to pg/ml (Original: WT = pg/ml; OB = ng/ml)
Metabolome_blood(41:80,1) = Metabolome_blood(41:80,1)*1000;

%% Parameters
num_Genotype = 2;
index_WT = 1;
index_OB = 2;
num_TissueType = 1;
index_blood = 1;
num_TimePoint = 8;
TimePoint = [0,2,4,6,8,12,16,24];
col_num_blood = [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5];
num_Replicate = 5;
Group = {'WT0','WT0','WT0','WT0','WT0','WT2','WT2','WT2','WT2','WT2','WT4','WT4','WT4','WT4','WT4','WT6','WT6','WT6','WT6','WT6','WT8','WT8','WT8','WT8','WT8','WT12','WT12','WT12','WT12','WT12','WT16','WT16','WT16','WT16','WT16','WT24','WT24','WT24','WT24','WT24','OB0','OB0','OB0','OB0','OB0','OB2','OB2','OB2','OB2','OB2','OB4','OB4','OB4','OB4','OB4','OB6','OB6','OB6','OB6','OB6','OB8','OB8','OB8','OB8','OB8','OB12','OB12','OB12','OB12','OB12','OB16','OB16','OB16','OB16','OB16','OB24','OB24','OB24','OB24','OB24'};
threshold_q = 0.1;

%% Check Enough Data (Exist in >= half replicates)
is_EnoughData = false(num_TissueType,size(Info_MetabolomeName,1));
for i = 1:size(is_EnoughData,2)
    logic_EnoughData = CheckEnoughData(Metabolome_blood(:,i), col_num_blood);
    is_EnoughData(index_blood,i) = logic_EnoughData;
end

%%% Calculation
%%% Mean
mean_blood = NaN(length(col_num_blood),size(Info_MetabolomeName,1));
for i = 1:size(mean_blood,2)
    mean_blood(:,i) = GetMean(Metabolome_blood(:,i),col_num_blood);
end

%%% FC & log2FC (WT tX/t0, OB tX/t0, OB/WT tX)(x = 0,2,4,6,8,12,16,24)
FC_blood = [mean_blood(1:num_TimePoint,:)./mean_blood(1,:);mean_blood(num_TimePoint+1:num_TimePoint*2,:)./mean_blood(num_TimePoint+1,:);mean_blood(num_TimePoint+1:num_TimePoint*2,:)./mean_blood(1:num_TimePoint,:)];
log2FC_blood = log2(FC_blood);

%%% p values (One-way ANOVA)(WT;OB;OB/WT tX)
p_blood = NaN(num_Genotype,size(Info_MetabolomeName,1));
temp_num = sum(col_num_blood(1:num_TimePoint));
% Blood - WT
for i = 1:size(Info_MetabolomeName,1)
    p_blood(index_WT,i) = anova1(Metabolome_blood(1:temp_num,i),Group(1:temp_num),'off');
end
% Blood - OB
for i = 1:size(Info_MetabolomeName,1)
    p_blood(index_OB,i) = anova1(Metabolome_blood(temp_num+1:end,i),Group(temp_num+1:end),'off');
end
% Blood - OB/WT
for i = 1:size(Info_MetabolomeName,1)
    for j = 1:num_TimePoint
        if j == 1
            temp_list_WT = Metabolome_blood(1:col_num_blood(j),i);
            temp_list_OB = Metabolome_blood(temp_num+1:temp_num+col_num_blood(num_TimePoint+j),i);
        else
            temp_list_WT = Metabolome_blood(sum(col_num_blood(1:j-1))+1:sum(col_num_blood(1:j)),i);
            temp_list_OB = Metabolome_blood(sum(col_num_blood(1:num_TimePoint+j-1))+1:sum(col_num_blood(1:num_TimePoint+j)),i);
        end
        p_blood(num_Genotype+j,i) = anova1([temp_list_WT,temp_list_OB],{'WT','OB'},'off');
    end
end

%%% q values for those with enough data (WT t0-tX;OB t0-tX;OB/WT tX)(x = 0,2,4,6,8,12,16,24)
%%% ('Method','polynomial','Lambda',0.01:0.01:0.90)or('BHFDR',true)
q_blood = NaN(size(p_blood));
for i = 1:size(p_blood,1)
    q_blood(i,is_EnoughData) = mafdr(p_blood(i,is_EnoughData),'BHFDR',true);
end

%%% AUC (log2(data/geometric mean of WT_0 and ob_0))
geomean_time0 = geomean([mean_blood(1,:);mean_blood(num_TimePoint+1,:)],1);
log2normalize_bloodWT = log2(mean_blood(1:num_TimePoint,:)./geomean_time0);
AUC_curve_bloodWT = log2normalize_bloodWT - log2normalize_bloodWT(1,:);
log2normalize_bloodOB = log2(mean_blood(num_TimePoint+1:end,:)./geomean_time0);
AUC_curve_bloodOB = log2normalize_bloodOB - log2normalize_bloodOB(1,:);
AUC_bloodWT = trapz(TimePoint,AUC_curve_bloodWT,1);
AUC_bloodOB = trapz(TimePoint,AUC_curve_bloodOB,1);

%%% Responsive & Different (with Enough Data)
is_Responsive(index_WT,:) = (q_blood(index_WT,:)<threshold_q).*(AUC_bloodWT>0)*1 + (q_blood(index_WT,:)<threshold_q).*(AUC_bloodWT<0)*(-1);
is_Responsive(index_OB,:) = (q_blood(index_OB,:)<threshold_q).*(AUC_bloodOB>0)*1 + (q_blood(index_OB,:)<threshold_q).*(AUC_bloodOB<0)*(-1);
is_Different(1:num_TimePoint,:) = (q_blood(num_Genotype+1:end,:)<threshold_q).*(FC_blood(num_TimePoint*2+1:num_TimePoint*3,:)>1)*1 + (q_blood(num_Genotype+1:end,:)<threshold_q).*(FC_blood(num_TimePoint*2+1:num_TimePoint*3,:)<1)*(-1);
% Mixed as 2
is_TimeCourseDifferent = false(1,size(is_Different,2));
is_TimeCourseDifferent = is_TimeCourseDifferent + (sum(logical(is_Different))~=0).*(logical(~sum(is_Different==-1))).*(1);
is_TimeCourseDifferent = is_TimeCourseDifferent + (sum(logical(is_Different))~=0).*(logical(~sum(is_Different==1))).*(-1);
is_TimeCourseDifferent = is_TimeCourseDifferent + (sum(logical(is_Different))~=0).*(logical(sum(is_Different==1)).*logical(sum(is_Different==-1))).*(2);

%%% Calculate T half
% 1=fast,2=medium,3=slow
% {"<8 h","8-16 h",">16 h"}
THalf_blood = NaN(num_Genotype,size(Info_MetabolomeName,1));
for i = 1:size(THalf_blood,1)
    for j = 1:size(THalf_blood,2)
        temp_y = mean_blood((i-1)*num_TimePoint+1:num_TimePoint*i,j);
        if is_Responsive(i,j) == 1
            temp_ymax = max(temp_y);
            temp_y0 = mean([temp_y(1),temp_ymax]);
            temp_index_max = find(temp_y == temp_ymax);
            for k = 1:temp_index_max-1
                if (temp_y(k)<=temp_y0)&&(temp_y0<=temp_y(k+1))
                    temp_THalf = GetTHalf(TimePoint(k),temp_y(k),TimePoint(k+1),temp_y(k+1),temp_y0);
                    THalf_blood(i,j) = temp_THalf;
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
                    THalf_blood(i,j) = temp_THalf;
                    break
                end
            end
        end
    end
end
THalfType_blood = (THalf_blood<8)*1+(THalf_blood>=8).*(THalf_blood<=16)*2+(THalf_blood>16)*3;

%% Figure
%%% Insulin
Insulin_WT_matrix = reshape(Metabolome_blood(1:40,1),num_Replicate,num_TimePoint);
Insulin_OB_matrix = reshape(Metabolome_blood(41:80,1),num_Replicate,num_TimePoint);
Insulin_WT_mean = mean(Insulin_WT_matrix,1,'omitnan');
Insulin_OB_mean = mean(Insulin_OB_matrix,1,'omitnan');
Insulin_WT_error = std(Insulin_WT_matrix,0,1,'omitnan');
Insulin_OB_error = std(Insulin_OB_matrix,0,1,'omitnan');

f1 = figure(1);
hold on;
for i = 1:size(Insulin_WT_matrix,1)
    scatter(TimePoint,Insulin_WT_matrix(i,:),[],[0.8 0.8 1],'filled');
end
for i = 1:size(Insulin_OB_matrix,1)
    scatter(TimePoint,Insulin_OB_matrix(i,:),[],[1 0.8 0.8],'filled');
end
p1 = errorbar(TimePoint,Insulin_WT_mean,Insulin_WT_error,'Color',[0 0 1],'LineWidth',1);
p2 = errorbar(TimePoint,Insulin_OB_mean,Insulin_OB_error,'Color',[1 0 0],'LineWidth',1);
%text(8,0.93,'*','FontSize',20,'Color','k','HorizontalAlignment','center');          % Manual

xlim([0,24]);
xticks(TimePoint);
xlabel("Time Course (h)");
%ylim([0.8,1]);                                                                      % Manual
ylabel("Blood Insulin (pg/ml)");
legend([p1,p2],["WT","\it{ob}\rm{/}\it{ob}"]);
f1.Position = [100 500 220 200];

title(['WT: p = ',num2str(p_blood(1,1)),sprintf('\n'),'\it{ob}\rm\bf{/}\it{ob}\rm\bf{: p = }',num2str(p_blood(2,1))])

%%
exportgraphics(f1,['./Fig/','Blood Insulin','.pdf'],'ContentType','vector');

%% Figure
%%% Glucose
Glucose_WT_matrix = reshape(Metabolome_blood(1:40,2),num_Replicate,num_TimePoint);
Glucose_OB_matrix = reshape(Metabolome_blood(41:80,2),num_Replicate,num_TimePoint);
Glucose_WT_mean = mean(Glucose_WT_matrix,1,'omitnan');
Glucose_OB_mean = mean(Glucose_OB_matrix,1,'omitnan');
Glucose_WT_error = std(Glucose_WT_matrix,0,1,'omitnan');
Glucose_OB_error = std(Glucose_OB_matrix,0,1,'omitnan');

f2 = figure(2);
hold on;
for i = 1:size(Glucose_WT_matrix,1)
    scatter(TimePoint,Glucose_WT_matrix(i,:),[],[0.8 0.8 1],'filled');
end
for i = 1:size(Glucose_OB_matrix,1)
    scatter(TimePoint,Glucose_OB_matrix(i,:),[],[1 0.8 0.8],'filled');
end
p1 = errorbar(TimePoint,Glucose_WT_mean,Glucose_WT_error,'Color',[0 0 1],'LineWidth',1);
p2 = errorbar(TimePoint,Glucose_OB_mean,Glucose_OB_error,'Color',[1 0 0],'LineWidth',1);
%text(8,0.93,'*','FontSize',20,'Color','k','HorizontalAlignment','center');          % Manual

xlim([0,24]);
xticks(TimePoint);
xlabel("Time Course (h)");
%ylim([0.8,1]);                                                                      % Manual
ylabel("Blood Glucose (mg/dl)");
legend([p1,p2],["WT","\it{ob}\rm{/}\it{ob}"]);
f2.Position = [100 500 220 200];

title(['WT: p = ',num2str(p_blood(1,2)),sprintf('\n'),'\it{ob}\rm\bf{/}\it{ob}\rm\bf{: p = }',num2str(p_blood(2,2))])

%%
exportgraphics(f2,['./Fig/','Blood Glucose','.pdf'],'ContentType','vector');

%% %%%%%%%%%%%%%%%%%%%%%%%%%

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
