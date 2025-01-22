%% Load files
load('./Resource/WB_StarvationLiver.mat');

%%
Info_ProteinName_percentage = cellfun(@(s) ['% ' s], Info_ProteinName(mod([1:length(Info_ProteinName)],2) == 0), 'UniformOutput', false);

Data_Total = Data(:,mod([1:length(Info_ProteinName)],2) == 1);
Data_Phospho = Data(:,mod([1:length(Info_ProteinName)],2) == 0);
Data_Percentage = Data_Phospho./Data_Total;

%% Parameters
num_Genotype = 2;
index_WT = 1;
index_OB = 2;
num_TissueType = 1;
index_liver = 1;
num_TimePoint = 8;
TimePoint = [0,2,4,6,8,12,16,24];
col_num_liver = [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5];
Group = {'WT0','WT0','WT0','WT0','WT0','WT2','WT2','WT2','WT2','WT2','WT4','WT4','WT4','WT4','WT4','WT6','WT6','WT6','WT6','WT6','WT8','WT8','WT8','WT8','WT8','WT12','WT12','WT12','WT12','WT12','WT16','WT16','WT16','WT16','WT16','WT24','WT24','WT24','WT24','WT24','OB0','OB0','OB0','OB0','OB0','OB2','OB2','OB2','OB2','OB2','OB4','OB4','OB4','OB4','OB4','OB6','OB6','OB6','OB6','OB6','OB8','OB8','OB8','OB8','OB8','OB12','OB12','OB12','OB12','OB12','OB16','OB16','OB16','OB16','OB16','OB24','OB24','OB24','OB24','OB24'};
threshold_q = 0.1;

%% Mean
mean_liver = NaN(length(col_num_liver),size(Data,2));
for i = 1:length(col_num_liver)
    mean_liver(i,:) = mean(Data([1:5]+5*(i-1),:),1);
end

mean_percentage = NaN(length(col_num_liver),size(Data_Percentage,2));
for i = 1:length(col_num_liver)
    mean_percentage(i,:) = mean(Data_Percentage([1:5]+5*(i-1),:),1);
end

%% p values (One-way ANOVA)(WT t0-tX;OB t0-tX;OB/WT tX)
p_liver = NaN(num_Genotype,size(Data,2));
temp_num = sum(col_num_liver(1:num_TimePoint));

%%% For total and phospho amount
% liver - WT
for i = 1:size(Data,2)
    p_liver(index_WT,i) = anova1(Data(1:temp_num,i),Group(1:temp_num),'off');
end
% liver - OB
for i = 1:size(Data,2)
    p_liver(index_OB,i) = anova1(Data(temp_num+1:end,i),Group(temp_num+1:end),'off');
end
% liver - OB/WT
for i = 1:size(Data,2)
    for j = 1:num_TimePoint
        if j == 1
            temp_list_WT = Data(1:col_num_liver(j),i);
            temp_list_OB = Data(temp_num+1:temp_num+col_num_liver(num_TimePoint+j),i);
        else
            temp_list_WT = Data(sum(col_num_liver(1:j-1))+1:sum(col_num_liver(1:j)),i);
            temp_list_OB = Data(sum(col_num_liver(1:num_TimePoint+j-1))+1:sum(col_num_liver(1:num_TimePoint+j)),i);
        end
        p_liver(num_Genotype+j,i) = anova1([temp_list_WT,temp_list_OB],{'WT','OB'},'off');
    end
end

%%% For percentage
% liver - WT
for i = 1:size(Data_Percentage,2)
    p_percentage(index_WT,i) = anova1(Data_Percentage(1:temp_num,i),Group(1:temp_num),'off');
end
% liver - OB
for i = 1:size(Data_Percentage,2)
    p_percentage(index_OB,i) = anova1(Data_Percentage(temp_num+1:end,i),Group(temp_num+1:end),'off');
end
% liver - OB/WT
for i = 1:size(Data_Percentage,2)
    for j = 1:num_TimePoint
        if j == 1
            temp_list_WT = Data_Percentage(1:col_num_liver(j),i);
            temp_list_OB = Data_Percentage(temp_num+1:temp_num+col_num_liver(num_TimePoint+j),i);
        else
            temp_list_WT = Data_Percentage(sum(col_num_liver(1:j-1))+1:sum(col_num_liver(1:j)),i);
            temp_list_OB = Data_Percentage(sum(col_num_liver(1:num_TimePoint+j-1))+1:sum(col_num_liver(1:num_TimePoint+j)),i);
        end
        p_percentage(num_Genotype+j,i) = anova1([temp_list_WT,temp_list_OB],{'WT','OB'},'off');
    end
end

%% q values (WT t0-tX;OB t0-tX;OB/WT tX)
%%% ('Method','polynomial','Lambda',0.01:0.01:0.90)or('BHFDR',true)

%%% For total and phospho amount
q_liver = NaN(size(p_liver));
for i = 1:size(p_liver,1)
    q_liver(i,:) = mafdr(p_liver(i,:),'BHFDR',true);
end
%%% For percentage
q_percentage = NaN(size(p_percentage));
for i = 1:size(p_percentage,1)
    q_percentage(i,:) = mafdr(p_percentage(i,:),'BHFDR',true);
end

%% AUC (log2(data/geometric mean of WT_0 and ob_0))
%%% For total and phospho amount
geomean_time0 = geomean([mean_liver(1,:);mean_liver(num_TimePoint+1,:)],1);
log2normalize_liverWT = log2(mean_liver(1:num_TimePoint,:)./geomean_time0);
AUC_curve_liverWT = log2normalize_liverWT - log2normalize_liverWT(1,:);
log2normalize_liverOB = log2(mean_liver(num_TimePoint+1:end,:)./geomean_time0);
AUC_curve_liverOB = log2normalize_liverOB - log2normalize_liverOB(1,:);
AUC_liverWT = trapz(TimePoint,AUC_curve_liverWT,1);
AUC_liverOB = trapz(TimePoint,AUC_curve_liverOB,1);

%%% For percentage
geomean_time0_percentage = geomean([mean_percentage(1,:);mean_percentage(num_TimePoint+1,:)],1);
log2normalize_liverWT_percentage = log2(mean_percentage(1:num_TimePoint,:)./geomean_time0_percentage);
AUC_curve_liverWT_percentage = log2normalize_liverWT_percentage - log2normalize_liverWT_percentage(1,:);
log2normalize_liverOB_percentage = log2(mean_percentage(num_TimePoint+1:end,:)./geomean_time0_percentage);
AUC_curve_liverOB_percentage = log2normalize_liverOB_percentage - log2normalize_liverOB_percentage(1,:);
AUC_liverWT_percentage = trapz(TimePoint,AUC_curve_liverWT_percentage,1);
AUC_liverOB_percentage = trapz(TimePoint,AUC_curve_liverOB_percentage,1);

%% Responsive & Different
%%% For total and phospho amount
is_Responsive(index_WT,:) = (q_liver(index_WT,:)<threshold_q).*(AUC_liverWT>0)*1 + (q_liver(index_WT,:)<threshold_q).*(AUC_liverWT<0)*(-1);
is_Responsive(index_OB,:) = (q_liver(index_OB,:)<threshold_q).*(AUC_liverOB>0)*1 + (q_liver(index_OB,:)<threshold_q).*(AUC_liverOB<0)*(-1);
is_Different(1:num_TimePoint,:) = (q_liver(num_Genotype+1:end,:)<threshold_q).*(mean_liver(1:num_TimePoint,:) < mean_liver(num_TimePoint+1:num_TimePoint*2,:))*1 + (q_liver(num_Genotype+1:end,:)<threshold_q).*(mean_liver(1:num_TimePoint,:) > mean_liver(num_TimePoint+1:num_TimePoint*2,:))*(-1);

%%% For percentage
is_Responsive_percentage(index_WT,:) = (q_percentage(index_WT,:)<threshold_q).*(AUC_liverWT_percentage>0)*1 + (q_percentage(index_WT,:)<threshold_q).*(AUC_liverWT_percentage<0)*(-1);
is_Responsive_percentage(index_OB,:) = (q_percentage(index_OB,:)<threshold_q).*(AUC_liverOB_percentage>0)*1 + (q_percentage(index_OB,:)<threshold_q).*(AUC_liverOB_percentage<0)*(-1);
is_Different_percentage(1:num_TimePoint,:) = (q_percentage(num_Genotype+1:end,:)<threshold_q).*(mean_percentage(1:num_TimePoint,:) < mean_percentage(num_TimePoint+1:num_TimePoint*2,:))*1 + (q_percentage(num_Genotype+1:end,:)<threshold_q).*(mean_percentage(1:num_TimePoint,:) > mean_percentage(num_TimePoint+1:num_TimePoint*2,:))*(-1);

%% Calculate T half
% 1=fast,2=medium,3=slow
% {"<8 h","8-16 h",">16 h"}
THalf_liver = NaN(num_Genotype,size(Info_ProteinName,1));
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

THalf_liver_percentage = NaN(num_Genotype,size(Info_ProteinName_percentage,1));
for i = 1:size(THalf_liver_percentage,1)
    for j = 1:size(THalf_liver_percentage,2)
        temp_y = mean_percentage((i-1)*num_TimePoint+1:num_TimePoint*i,j);
        if is_Responsive_percentage(i,j) == 1
            temp_ymax = max(temp_y);
            temp_y0 = mean([temp_y(1),temp_ymax]);
            temp_index_max = find(temp_y == temp_ymax);
            for k = 1:temp_index_max-1
                if (temp_y(k)<=temp_y0)&&(temp_y0<=temp_y(k+1))
                    temp_THalf = GetTHalf(TimePoint(k),temp_y(k),TimePoint(k+1),temp_y(k+1),temp_y0);
                    THalf_liver_percentage(i,j) = temp_THalf;
                    break
                end
            end
        elseif is_Responsive_percentage(i,j) == -1
            temp_ymin = min(temp_y);
            temp_y0 = mean([temp_y(1),temp_ymin]);
            temp_index_min = find(temp_y == temp_ymin);
            for k = 1:temp_index_min-1
                if (temp_y(k)>=temp_y0)&&(temp_y0>=temp_y(k+1))
                    temp_THalf = GetTHalf(TimePoint(k),temp_y(k),TimePoint(k+1),temp_y(k+1),temp_y0);
                    THalf_liver_percentage(i,j) = temp_THalf;
                    break
                end
            end
        end
    end
end
THalfType_liver_percentage = (THalf_liver_percentage<8)*1+(THalf_liver_percentage>=8).*(THalf_liver_percentage<=16)*2+(THalf_liver_percentage>16)*3;

%% %%%%%%%%%%%%%%%%%%%%%%%%%
%% Function
function THalf = GetTHalf(x1,y1,x2,y2,y0)
if y0 == y1
    THalf = x1;
elseif y0 == y2
    THalf = x2;
else
    THalf = x2 + (x2-x1)*(y0-y2)/(y2-y1);
end
end
