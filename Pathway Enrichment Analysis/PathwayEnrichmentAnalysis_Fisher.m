%% Load files
load('./Resource/PathwayEnrichmentAnalysis-Responsive&Dominant.mat');
% AllEnough 2 AllRes 2 PathEnough 2 PathRes 2 AllPath

%% Parameter
q_threshold = 0.01;

%% Responsive
%%% Calculate p value
p_WT_increase_T = NaN(length(Info_PathwayList),1);
p_WT_decrease_T = NaN(length(Info_PathwayList),1);
p_OB_increase_T = NaN(length(Info_PathwayList),1);
p_OB_decrease_T = NaN(length(Info_PathwayList),1);
p_WT_increase_P = NaN(length(Info_PathwayList),1);
p_WT_decrease_P = NaN(length(Info_PathwayList),1);
p_OB_increase_P = NaN(length(Info_PathwayList),1);
p_OB_decrease_P = NaN(length(Info_PathwayList),1);

for i = 1:length(Info_PathwayList)
    [~,p_WT_increase_T(i,1),~] = fishertest([Input_WT_increase_T(i,4),Input_WT_increase_T(i,3)-Input_WT_increase_T(i,4);Input_WT_increase_T(i,2)-Input_WT_increase_T(i,4),Input_WT_increase_T(i,1)-Input_WT_increase_T(i,2)-Input_WT_increase_T(i,3)+Input_WT_increase_T(i,4)]);
    if Input_WT_increase_T(i,4) == 0
        p_WT_increase_T(i,1) = 1;
    end
end
for i = 1:length(Info_PathwayList)
    [~,p_WT_decrease_T(i,1),~] = fishertest([Input_WT_decrease_T(i,4),Input_WT_decrease_T(i,3)-Input_WT_decrease_T(i,4);Input_WT_decrease_T(i,2)-Input_WT_decrease_T(i,4),Input_WT_decrease_T(i,1)-Input_WT_decrease_T(i,2)-Input_WT_decrease_T(i,3)+Input_WT_decrease_T(i,4)]);
    if Input_WT_decrease_T(i,4) == 0
        p_WT_decrease_T(i,1) = 1;
    end
end
for i = 1:length(Info_PathwayList)
    [~,p_OB_increase_T(i,1),~] = fishertest([Input_OB_increase_T(i,4),Input_OB_increase_T(i,3)-Input_OB_increase_T(i,4);Input_OB_increase_T(i,2)-Input_OB_increase_T(i,4),Input_OB_increase_T(i,1)-Input_OB_increase_T(i,2)-Input_OB_increase_T(i,3)+Input_OB_increase_T(i,4)]);
    if Input_OB_increase_T(i,4) == 0
        p_OB_increase_T(i,1) = 1;
    end
end
for i = 1:length(Info_PathwayList)
    [~,p_OB_decrease_T(i,1),~] = fishertest([Input_OB_decrease_T(i,4),Input_OB_decrease_T(i,3)-Input_OB_decrease_T(i,4);Input_OB_decrease_T(i,2)-Input_OB_decrease_T(i,4),Input_OB_decrease_T(i,1)-Input_OB_decrease_T(i,2)-Input_OB_decrease_T(i,3)+Input_OB_decrease_T(i,4)]);
    if Input_OB_decrease_T(i,4) == 0
        p_OB_decrease_T(i,1) = 1;
    end
end
for i = 1:length(Info_PathwayList)
    [~,p_WT_increase_P(i,1),~] = fishertest([Input_WT_increase_P(i,4),Input_WT_increase_P(i,3)-Input_WT_increase_P(i,4);Input_WT_increase_P(i,2)-Input_WT_increase_P(i,4),Input_WT_increase_P(i,1)-Input_WT_increase_P(i,2)-Input_WT_increase_P(i,3)+Input_WT_increase_P(i,4)]);
    if Input_WT_increase_P(i,4) == 0
        p_WT_increase_P(i,1) = 1;
    end
end
for i = 1:length(Info_PathwayList)
    [~,p_WT_decrease_P(i,1),~] = fishertest([Input_WT_decrease_P(i,4),Input_WT_decrease_P(i,3)-Input_WT_decrease_P(i,4);Input_WT_decrease_P(i,2)-Input_WT_decrease_P(i,4),Input_WT_decrease_P(i,1)-Input_WT_decrease_P(i,2)-Input_WT_decrease_P(i,3)+Input_WT_decrease_P(i,4)]);
    if Input_WT_decrease_P(i,4) == 0
        p_WT_decrease_P(i,1) = 1;
    end
end
for i = 1:length(Info_PathwayList)
    [~,p_OB_increase_P(i,1),~] = fishertest([Input_OB_increase_P(i,4),Input_OB_increase_P(i,3)-Input_OB_increase_P(i,4);Input_OB_increase_P(i,2)-Input_OB_increase_P(i,4),Input_OB_increase_P(i,1)-Input_OB_increase_P(i,2)-Input_OB_increase_P(i,3)+Input_OB_increase_P(i,4)]);
    if Input_OB_increase_P(i,4) == 0
        p_OB_increase_P(i,1) = 1;
    end
end
for i = 1:length(Info_PathwayList)
    [~,p_OB_decrease_P(i,1),~] = fishertest([Input_OB_decrease_P(i,4),Input_OB_decrease_P(i,3)-Input_OB_decrease_P(i,4);Input_OB_decrease_P(i,2)-Input_OB_decrease_P(i,4),Input_OB_decrease_P(i,1)-Input_OB_decrease_P(i,2)-Input_OB_decrease_P(i,3)+Input_OB_decrease_P(i,4)]);
    if Input_OB_decrease_P(i,4) == 0
        p_OB_decrease_P(i,1) = 1;
    end
end

%%% Calculate q value
q_WT_increase_T = mafdr(p_WT_increase_T,'BHFDR',true);
q_WT_decrease_T = mafdr(p_WT_decrease_T,'BHFDR',true);
q_OB_increase_T = mafdr(p_OB_increase_T,'BHFDR',true);
q_OB_decrease_T = mafdr(p_OB_decrease_T,'BHFDR',true);
q_WT_increase_P = mafdr(p_WT_increase_P,'BHFDR',true);
q_WT_decrease_P = mafdr(p_WT_decrease_P,'BHFDR',true);
q_OB_increase_P = mafdr(p_OB_increase_P,'BHFDR',true);
q_OB_decrease_P = mafdr(p_OB_decrease_P,'BHFDR',true);

%%% Result
Result_WT_increase_T = Info_PathwayList(q_WT_increase_T < q_threshold,:);
Result_WT_decrease_T = Info_PathwayList(q_WT_decrease_T < q_threshold,:);
Result_OB_increase_T = Info_PathwayList(q_OB_increase_T < q_threshold,:);
Result_OB_decrease_T = Info_PathwayList(q_OB_decrease_T < q_threshold,:);
Result_WT_increase_P = Info_PathwayList(q_WT_increase_P < q_threshold,:);
Result_WT_decrease_P = Info_PathwayList(q_WT_decrease_P < q_threshold,:);
Result_OB_increase_P = Info_PathwayList(q_OB_increase_P < q_threshold,:);
Result_OB_decrease_P = Info_PathwayList(q_OB_decrease_P < q_threshold,:);

%%
%{
%% THalf
%%% Calculate p value
p_WT_1_T = NaN(length(Info_PathwayList),1);
p_WT_2_T = NaN(length(Info_PathwayList),1);
p_WT_3_T = NaN(length(Info_PathwayList),1);
p_OB_1_T = NaN(length(Info_PathwayList),1);
p_OB_2_T = NaN(length(Info_PathwayList),1);
p_OB_3_T = NaN(length(Info_PathwayList),1);
p_WT_1_P = NaN(length(Info_PathwayList),1);
p_WT_2_P = NaN(length(Info_PathwayList),1);
p_WT_3_P = NaN(length(Info_PathwayList),1);
p_OB_1_P = NaN(length(Info_PathwayList),1);
p_OB_2_P = NaN(length(Info_PathwayList),1);
p_OB_3_P = NaN(length(Info_PathwayList),1);

for i = 1:length(Info_PathwayList)
    [~,p_WT_1_T(i,1),~] = fishertest([Input_WT_1_T(i,4),Input_WT_1_T(i,3)-Input_WT_1_T(i,4);Input_WT_1_T(i,2)-Input_WT_1_T(i,4),Input_WT_1_T(i,1)-Input_WT_1_T(i,2)-Input_WT_1_T(i,3)+Input_WT_1_T(i,4)]);
    if Input_WT_1_T(i,4) == 0
        p_WT_1_T(i,1) = 1;
    end
end
for i = 1:length(Info_PathwayList)
    [~,p_WT_2_T(i,1),~] = fishertest([Input_WT_2_T(i,4),Input_WT_2_T(i,3)-Input_WT_2_T(i,4);Input_WT_2_T(i,2)-Input_WT_2_T(i,4),Input_WT_2_T(i,1)-Input_WT_2_T(i,2)-Input_WT_2_T(i,3)+Input_WT_2_T(i,4)]);
    if Input_WT_2_T(i,4) == 0
        p_WT_2_T(i,1) = 1;
    end
end
for i = 1:length(Info_PathwayList)
    [~,p_WT_3_T(i,1),~] = fishertest([Input_WT_3_T(i,4),Input_WT_3_T(i,3)-Input_WT_3_T(i,4);Input_WT_3_T(i,2)-Input_WT_3_T(i,4),Input_WT_3_T(i,1)-Input_WT_3_T(i,2)-Input_WT_3_T(i,3)+Input_WT_3_T(i,4)]);
    if Input_WT_3_T(i,4) == 0
        p_WT_3_T(i,1) = 1;
    end
end
for i = 1:length(Info_PathwayList)
    [~,p_OB_1_T(i,1),~] = fishertest([Input_OB_1_T(i,4),Input_OB_1_T(i,3)-Input_OB_1_T(i,4);Input_OB_1_T(i,2)-Input_OB_1_T(i,4),Input_OB_1_T(i,1)-Input_OB_1_T(i,2)-Input_OB_1_T(i,3)+Input_OB_1_T(i,4)]);
    if Input_OB_1_T(i,4) == 0
        p_OB_1_T(i,1) = 1;
    end
end
for i = 1:length(Info_PathwayList)
    [~,p_OB_2_T(i,1),~] = fishertest([Input_OB_2_T(i,4),Input_OB_2_T(i,3)-Input_OB_2_T(i,4);Input_OB_2_T(i,2)-Input_OB_2_T(i,4),Input_OB_2_T(i,1)-Input_OB_2_T(i,2)-Input_OB_2_T(i,3)+Input_OB_2_T(i,4)]);
    if Input_OB_2_T(i,4) == 0
        p_OB_2_T(i,1) = 1;
    end
end
for i = 1:length(Info_PathwayList)
    [~,p_OB_3_T(i,1),~] = fishertest([Input_OB_3_T(i,4),Input_OB_3_T(i,3)-Input_OB_3_T(i,4);Input_OB_3_T(i,2)-Input_OB_3_T(i,4),Input_OB_3_T(i,1)-Input_OB_3_T(i,2)-Input_OB_3_T(i,3)+Input_OB_3_T(i,4)]);
    if Input_OB_3_T(i,4) == 0
        p_OB_3_T(i,1) = 1;
    end
end
for i = 1:length(Info_PathwayList)
    [~,p_WT_1_P(i,1),~] = fishertest([Input_WT_1_P(i,4),Input_WT_1_P(i,3)-Input_WT_1_P(i,4);Input_WT_1_P(i,2)-Input_WT_1_P(i,4),Input_WT_1_P(i,1)-Input_WT_1_P(i,2)-Input_WT_1_P(i,3)+Input_WT_1_P(i,4)]);
    if Input_WT_1_P(i,4) == 0
        p_WT_1_P(i,1) = 1;
    end
end
for i = 1:length(Info_PathwayList)
    [~,p_WT_2_P(i,1),~] = fishertest([Input_WT_2_P(i,4),Input_WT_2_P(i,3)-Input_WT_2_P(i,4);Input_WT_2_P(i,2)-Input_WT_2_P(i,4),Input_WT_2_P(i,1)-Input_WT_2_P(i,2)-Input_WT_2_P(i,3)+Input_WT_2_P(i,4)]);
    if Input_WT_2_P(i,4) == 0
        p_WT_2_P(i,1) = 1;
    end
end
for i = 1:length(Info_PathwayList)
    [~,p_WT_3_P(i,1),~] = fishertest([Input_WT_3_P(i,4),Input_WT_3_P(i,3)-Input_WT_3_P(i,4);Input_WT_3_P(i,2)-Input_WT_3_P(i,4),Input_WT_3_P(i,1)-Input_WT_3_P(i,2)-Input_WT_3_P(i,3)+Input_WT_3_P(i,4)]);
    if Input_WT_3_P(i,4) == 0
        p_WT_3_P(i,1) = 1;
    end
end
for i = 1:length(Info_PathwayList)
    [~,p_OB_1_P(i,1),~] = fishertest([Input_OB_1_P(i,4),Input_OB_1_P(i,3)-Input_OB_1_P(i,4);Input_OB_1_P(i,2)-Input_OB_1_P(i,4),Input_OB_1_P(i,1)-Input_OB_1_P(i,2)-Input_OB_1_P(i,3)+Input_OB_1_P(i,4)]);
    if Input_OB_1_P(i,4) == 0
        p_OB_1_P(i,1) = 1;
    end
end
for i = 1:length(Info_PathwayList)
    [~,p_OB_2_P(i,1),~] = fishertest([Input_OB_2_P(i,4),Input_OB_2_P(i,3)-Input_OB_2_P(i,4);Input_OB_2_P(i,2)-Input_OB_2_P(i,4),Input_OB_2_P(i,1)-Input_OB_2_P(i,2)-Input_OB_2_P(i,3)+Input_OB_2_P(i,4)]);
    if Input_OB_2_P(i,4) == 0
        p_OB_2_P(i,1) = 1;
    end
end
for i = 1:length(Info_PathwayList)
    [~,p_OB_3_P(i,1),~] = fishertest([Input_OB_3_P(i,4),Input_OB_3_P(i,3)-Input_OB_3_P(i,4);Input_OB_3_P(i,2)-Input_OB_3_P(i,4),Input_OB_3_P(i,1)-Input_OB_3_P(i,2)-Input_OB_3_P(i,3)+Input_OB_3_P(i,4)]);
    if Input_OB_3_P(i,4) == 0
        p_OB_3_P(i,1) = 1;
    end
end

%%% Calculate q value
q_WT_1_T = mafdr(p_WT_1_T,'BHFDR',true);
q_WT_2_T = mafdr(p_WT_2_T,'BHFDR',true);
q_WT_3_T = mafdr(p_WT_3_T,'BHFDR',true);
q_OB_1_T = mafdr(p_OB_1_T,'BHFDR',true);
q_OB_2_T = mafdr(p_OB_2_T,'BHFDR',true);
q_OB_3_T = mafdr(p_OB_3_T,'BHFDR',true);
q_WT_1_P = mafdr(p_WT_1_P,'BHFDR',true);
q_WT_2_P = mafdr(p_WT_2_P,'BHFDR',true);
q_WT_3_P = mafdr(p_WT_3_P,'BHFDR',true);
q_OB_1_P = mafdr(p_OB_1_P,'BHFDR',true);
q_OB_2_P = mafdr(p_OB_2_P,'BHFDR',true);
q_OB_3_P = mafdr(p_OB_3_P,'BHFDR',true);

%%% Result
Result_WT_1_T = Info_PathwayList(q_WT_1_T < q_threshold,:);
Result_WT_2_T = Info_PathwayList(q_WT_2_T < q_threshold,:);
Result_WT_3_T = Info_PathwayList(q_WT_3_T < q_threshold,:);
Result_OB_1_T = Info_PathwayList(q_OB_1_T < q_threshold,:);
Result_OB_2_T = Info_PathwayList(q_OB_2_T < q_threshold,:);
Result_OB_3_T = Info_PathwayList(q_OB_3_T < q_threshold,:);
Result_WT_1_P = Info_PathwayList(q_WT_1_P < q_threshold,:);
Result_WT_2_P = Info_PathwayList(q_WT_2_P < q_threshold,:);
Result_WT_3_P = Info_PathwayList(q_WT_3_P < q_threshold,:);
Result_OB_1_P = Info_PathwayList(q_OB_1_P < q_threshold,:);
Result_OB_2_P = Info_PathwayList(q_OB_2_P < q_threshold,:);
Result_OB_3_P = Info_PathwayList(q_OB_3_P < q_threshold,:);
%}

%% Dominant
%%% Calculate p value
p_WT_T = NaN(length(Info_PathwayList),1);
p_OB_T = NaN(length(Info_PathwayList),1);
p_mix_T = NaN(length(Info_PathwayList),1);
p_WT_P = NaN(length(Info_PathwayList),1);
p_OB_P = NaN(length(Info_PathwayList),1);
p_mix_P = NaN(length(Info_PathwayList),1);

for i = 1:length(Info_PathwayList)
    [~,p_WT_T(i,1),~] = fishertest([Input_Dominant_WT_T(i,4),Input_Dominant_WT_T(i,3)-Input_Dominant_WT_T(i,4);Input_Dominant_WT_T(i,2)-Input_Dominant_WT_T(i,4),Input_Dominant_WT_T(i,1)-Input_Dominant_WT_T(i,2)-Input_Dominant_WT_T(i,3)+Input_Dominant_WT_T(i,4)]);
    if Input_Dominant_WT_T(i,4) == 0
        p_WT_T(i,1) = 1;
    end
end
for i = 1:length(Info_PathwayList)
    [~,p_OB_T(i,1),~] = fishertest([Input_Dominant_OB_T(i,4),Input_Dominant_OB_T(i,3)-Input_Dominant_OB_T(i,4);Input_Dominant_OB_T(i,2)-Input_Dominant_OB_T(i,4),Input_Dominant_OB_T(i,1)-Input_Dominant_OB_T(i,2)-Input_Dominant_OB_T(i,3)+Input_Dominant_OB_T(i,4)]);
    if Input_Dominant_OB_T(i,4) == 0
        p_OB_T(i,1) = 1;
    end
end
for i = 1:length(Info_PathwayList)
    [~,p_mix_T(i,1),~] = fishertest([Input_Dominant_mix_T(i,4),Input_Dominant_mix_T(i,3)-Input_Dominant_mix_T(i,4);Input_Dominant_mix_T(i,2)-Input_Dominant_mix_T(i,4),Input_Dominant_mix_T(i,1)-Input_Dominant_mix_T(i,2)-Input_Dominant_mix_T(i,3)+Input_Dominant_mix_T(i,4)]);
    if Input_Dominant_mix_T(i,4) == 0
        p_mix_T(i,1) = 1;
    end
end
for i = 1:length(Info_PathwayList)
    [~,p_WT_P(i,1),~] = fishertest([Input_Dominant_WT_P(i,4),Input_Dominant_WT_P(i,3)-Input_Dominant_WT_P(i,4);Input_Dominant_WT_P(i,2)-Input_Dominant_WT_P(i,4),Input_Dominant_WT_P(i,1)-Input_Dominant_WT_P(i,2)-Input_Dominant_WT_P(i,3)+Input_Dominant_WT_P(i,4)]);
    if Input_Dominant_WT_P(i,4) == 0
        p_WT_P(i,1) = 1;
    end
end
for i = 1:length(Info_PathwayList)
    [~,p_OB_P(i,1),~] = fishertest([Input_Dominant_OB_P(i,4),Input_Dominant_OB_P(i,3)-Input_Dominant_OB_P(i,4);Input_Dominant_OB_P(i,2)-Input_Dominant_OB_P(i,4),Input_Dominant_OB_P(i,1)-Input_Dominant_OB_P(i,2)-Input_Dominant_OB_P(i,3)+Input_Dominant_OB_P(i,4)]);
    if Input_Dominant_OB_P(i,4) == 0
        p_OB_P(i,1) = 1;
    end
end
for i = 1:length(Info_PathwayList)
    [~,p_mix_P(i,1),~] = fishertest([Input_Dominant_mix_P(i,4),Input_Dominant_mix_P(i,3)-Input_Dominant_mix_P(i,4);Input_Dominant_mix_P(i,2)-Input_Dominant_mix_P(i,4),Input_Dominant_mix_P(i,1)-Input_Dominant_mix_P(i,2)-Input_Dominant_mix_P(i,3)+Input_Dominant_mix_P(i,4)]);
    if Input_Dominant_mix_P(i,4) == 0
        p_mix_P(i,1) = 1;
    end
end

%%% Calculate q value
q_WT_T = mafdr(p_WT_T,'BHFDR',true);
q_OB_T = mafdr(p_OB_T,'BHFDR',true);
q_mix_T = mafdr(p_mix_T,'BHFDR',true);
q_WT_P = mafdr(p_WT_P,'BHFDR',true);
q_OB_P = mafdr(p_OB_P,'BHFDR',true);
q_mix_P = mafdr(p_mix_P,'BHFDR',true);

%%% Result
Result_WT_T = Info_PathwayList(q_WT_T < q_threshold,:);
Result_OB_T = Info_PathwayList(q_OB_T < q_threshold,:);
Result_mix_T = Info_PathwayList(q_mix_T < q_threshold,:);
Result_WT_P = Info_PathwayList(q_WT_P < q_threshold,:);
Result_OB_P = Info_PathwayList(q_OB_P < q_threshold,:);
Result_mix_P = Info_PathwayList(q_mix_P < q_threshold,:);

%%