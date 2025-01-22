%% Load files
load('./Resource/WB_OGTTMuscle.mat');

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

%% Mean
mean_muscle = NaN(length(col_num_muscle),size(Data,2));
for i = 1:length(col_num_muscle)
    mean_muscle(i,:) = mean(Data([1:5]+5*(i-1),:),1);
end

%% p values (One-way ANOVA)(WT t0-tX;OB t0-tX;OB/WT tX)
p_muscle = NaN(num_Genotype,size(Data,2));
temp_num = sum(col_num_muscle(1:num_TimePoint));

%%% For total and phospho amount
% muscle - WT
for i = 1:size(Data,2)
    p_muscle(index_WT,i) = anova1(Data(1:temp_num,i),Group(1:temp_num),'off');
end
% muscle - OB
for i = 1:size(Data,2)
    p_muscle(index_OB,i) = anova1(Data(temp_num+1:end,i),Group(temp_num+1:end),'off');
end
% muscle - OB/WT
for i = 1:size(Data,2)
    for j = 1:num_TimePoint
        if j == 1
            temp_list_WT = Data(1:col_num_muscle(j),i);
            temp_list_OB = Data(temp_num+1:temp_num+col_num_muscle(num_TimePoint+j),i);
        else
            temp_list_WT = Data(sum(col_num_muscle(1:j-1))+1:sum(col_num_muscle(1:j)),i);
            temp_list_OB = Data(sum(col_num_muscle(1:num_TimePoint+j-1))+1:sum(col_num_muscle(1:num_TimePoint+j)),i);
        end
        p_muscle(num_Genotype+j,i) = anova1([temp_list_WT,temp_list_OB],{'WT','OB'},'off');
    end
end

%% q values (WT t0-tX;OB t0-tX;OB/WT tX)
%%% ('Method','polynomial','Lambda',0.01:0.01:0.90)or('BHFDR',true)

%%% For total and phospho amount
q_muscle = NaN(size(p_muscle));
for i = 1:size(p_muscle,1)
    q_muscle(i,:) = mafdr(p_muscle(i,:),'BHFDR',true);
end

%% AUC (log2(data/geometric mean of WT_0 and ob_0))
%%% For total and phospho amount
geomean_time0 = geomean([mean_muscle(1,:);mean_muscle(num_TimePoint+1,:)],1);
log2normalize_muscleWT = log2(mean_muscle(1:num_TimePoint,:)./geomean_time0);
AUC_curve_muscleWT = log2normalize_muscleWT - log2normalize_muscleWT(1,:);
log2normalize_muscleOB = log2(mean_muscle(num_TimePoint+1:end,:)./geomean_time0);
AUC_curve_muscleOB = log2normalize_muscleOB - log2normalize_muscleOB(1,:);
AUC_muscleWT = trapz(TimePoint,AUC_curve_muscleWT,1);
AUC_muscleOB = trapz(TimePoint,AUC_curve_muscleOB,1);

%% Responsive & Different
%%% For total and phospho amount
is_Responsive(index_WT,:) = (q_muscle(index_WT,:)<threshold_q).*(AUC_muscleWT>0)*1 + (q_muscle(index_WT,:)<threshold_q).*(AUC_muscleWT<0)*(-1);
is_Responsive(index_OB,:) = (q_muscle(index_OB,:)<threshold_q).*(AUC_muscleOB>0)*1 + (q_muscle(index_OB,:)<threshold_q).*(AUC_muscleOB<0)*(-1);
is_Different(1:num_TimePoint,:) = (q_muscle(num_Genotype+1:end,:)<threshold_q).*(mean_muscle(1:num_TimePoint,:) < mean_muscle(num_TimePoint+1:num_TimePoint*2,:))*1 + (q_muscle(num_Genotype+1:end,:)<threshold_q).*(mean_muscle(1:num_TimePoint,:) > mean_muscle(num_TimePoint+1:num_TimePoint*2,:))*(-1);

%% Mixed as 2
is_TimeCourseDifferent = false(1,size(is_Different,2));
is_TimeCourseDifferent = is_TimeCourseDifferent + (sum(logical(is_Different))~=0).*(logical(~sum(is_Different==-1))).*(1);
is_TimeCourseDifferent = is_TimeCourseDifferent + (sum(logical(is_Different))~=0).*(logical(~sum(is_Different==1))).*(-1);
is_TimeCourseDifferent = is_TimeCourseDifferent + (sum(logical(is_Different))~=0).*(logical(sum(is_Different==1)).*logical(sum(is_Different==-1))).*(2);

%% Calculate T half
% 1=fast,2=medium,3=slow
% {"<8 h","8-16 h",">16 h"}
THalf_muscle = NaN(num_Genotype,size(Info_ProteinName,1));
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
THalfType_muscle = (THalf_muscle<8)*1+(THalf_muscle>=8).*(THalf_muscle<=16)*2+(THalf_muscle>16)*3;

THalf_muscle_percentage = NaN(num_Genotype,size(Info_ProteinName_percentage,1));
for i = 1:size(THalf_muscle_percentage,1)
    for j = 1:size(THalf_muscle_percentage,2)
        temp_y = mean_percentage((i-1)*num_TimePoint+1:num_TimePoint*i,j);
        if is_Responsive_percentage(i,j) == 1
            temp_ymax = max(temp_y);
            temp_y0 = mean([temp_y(1),temp_ymax]);
            temp_index_max = find(temp_y == temp_ymax);
            for k = 1:temp_index_max-1
                if (temp_y(k)<=temp_y0)&&(temp_y0<=temp_y(k+1))
                    temp_THalf = GetTHalf(TimePoint(k),temp_y(k),TimePoint(k+1),temp_y(k+1),temp_y0);
                    THalf_muscle_percentage(i,j) = temp_THalf;
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
                    THalf_muscle_percentage(i,j) = temp_THalf;
                    break
                end
            end
        end
    end
end
THalfType_muscle_percentage = (THalf_muscle_percentage<8)*1+(THalf_muscle_percentage>=8).*(THalf_muscle_percentage<=16)*2+(THalf_muscle_percentage>16)*3;

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
