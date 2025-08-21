function ExtractData_StressMemoryProject
%function ExtractData_StressMemoryProject
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Project: Stress and memory-based decisions %
%%% (with Lars Schwabe) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Extracting data and arranging them for R %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Inputs (defaults):
%%% [none]
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% sebastian.gluth@uni-hamburg.de %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%specifying the folder in which the output is saved
data_folder = 'Data';

%check how many datasets there are
data_sets = dir(data_folder); %screen folder
%ensure that only files starting with "Vp" are included
delete_data_sets = [];
for x = 1:size(data_sets,1)
    if length(data_sets(x).name)<3
        delete_data_sets = [delete_data_sets,x];
    else
        if sum(ismember('Vp',data_sets(x).name(1:2)))<2
            delete_data_sets = [delete_data_sets,x];
        end
    end
end
data_sets(delete_data_sets) = []; %get rid of wrong files
nsubj = size(data_sets,1);

%create matrices for each sub-task that will be forwarded to R

%%% --- %%% --- %%% --- %%%
%%% --- First and Second Rating (3 rounds of all 60 stimuli = 180 trials)
%%% --- %%% --- %%% --- %%%
ntrials_R = 60;
data_R = zeros(nsubj*ntrials_R,13); %ATTENTION: information will be merged across runs
%the 9 variables for the RATING matrix are:
	%1. subject ID
    %2. image number
    %3. trial number in run 1
    %4. trial number in run 2
    %5. trial number in run 3
    %6. rating in run 1
    %7. rating in run 2
    %8. rating in run 3
    %9. average rating across runs 1 and 2
    %10. average rating across all runs
    %11. response time in run 1
    %12. response time in run 2
    %13. response time in run 3
    
%%% --- %%% --- %%% --- %%%
%%% --- Encoding (3 blocks * 4 repetitions * 16 items = 192 trials)
%%% --- %%% --- %%% --- %%%
ntrials_E = 48;
nitems_per_block = ntrials_E/3;
data_E = zeros(nsubj*ntrials_E,20); %ATTENTION: information will be merged across runs (+ recall)
%the 20 variables for the ENCODING matrix are:
	%1. subject ID
    %2. image number
    %3. block number (background color)
    %4. image grid number (ATTENTION: this refers to the grid position)
    %5. trial number in run 1 (encoding)
    %6. trial number in run 1 (encoding recall)
    %7. trial number in run 2 (encoding)
    %8. trial number in run 2 (encoding recall)
    %9. trial number in run 3 (encoding)
    %10. trial number in run 3 (encoding recall)
    %11. trial number in run 4 (encoding)
    %12. trial number in run 4 (encoding recall)
    %13. response time in recall run 1 (blank squares)
    %14. response time in recall run 1 (squares with images shown)
    %15. response time in recall run 2 (blank squares)
    %16. response time in recall run 2 (squares with images shown)
    %17. response time in recall run 3 (blank squares)
    %18. response time in recall run 3 (squares with images shown)
    %19. response time in recall run 4 (blank squares)
    %20. response time in recall run 4 (squares with images shown)
    
%%% --- %%% --- %%% --- %%%
%%% --- Decision (3 blocks * 40 decisions = 120 trials)
%%% --- %%% --- %%% --- %%%
ntrials_D = 120;
data_D = zeros(nsubj*ntrials_D,12);
%the 12 variables for the DECISION matrix are:
	%1. subject ID
    %2. block number (background color)
    %3. control block (=1) or not (=0) [images are shown in control block]
    %4. trial
    %5. left image number
    %6. right image number
    %7. average rating of left image from ratings 1 and 2
    %8. average rating of left image from ratings 1, 2, and 3
    %9. average rating of right image from ratings 1 and 2
    %10. average rating of right image from ratings 1, 2, and 3
    %11. choice (1 = left, 2 = right, 3 = invalid)
    %12. response time
    
%%% --- %%% --- %%% --- %%%
%%% --- Cued Recall (3 blocks * 16 items = 48 trials)
%%% --- %%% --- %%% --- %%%
ntrials_C = 48;
data_C = zeros(nsubj*ntrials_C,8);
%the 8 variables for the CUED_RECALL matrix are:
	%1. subject ID
    %2. block number (background color)
    %3. control block (=1) or not (=0) [images shown only during decision!]
    %4. trial
    %5. image number
    %6. response time (when saying the item's identity aloud)
    %7. confidence rating
    %8. response time of confidence rating

%go into single subjects to extract data
subject_ID = zeros(nsubj,1);
for s = 1:nsubj
    cdata = load([data_folder,'/',data_sets(s).name]); %read in .mat file
      
    %get subject ID from file name
    if isnan(str2double(data_sets(s).name(4)))
        subject_ID(s) = str2double(data_sets(s).name(3));
    elseif isnan(str2double(data_sets(s).name(5)))
        subject_ID(s) = str2double(data_sets(s).name(3:4));
    else
        subject_ID(s) = str2double(data_sets(s).name(3:5));
    end
    
    %%% --- %%% --- %%% --- %%%
    %%% --- First and Second Rating (3 rounds of all 60 stimuli = 180 trials)
    %%% --- %%% --- %%% --- %%%
    data_R_s = zeros(ntrials_R,size(data_R,2)); %only for current subject
    for t = 1:ntrials_R %loop over trials of 1st and 3rd run
        tdata = cdata.ImageRating1(t);
        data_R_s(t,[1:3,6,11]) = [subject_ID(s),tdata.imageNumber,t,tdata.response,tdata.onsetResponse-tdata.onsetDisplay];
    end
    for t = ntrials_R+1:ntrials_R*2 %loop over trials of 2nd run
        tdata = cdata.ImageRating1(t);
        %find trial in run 2 that matches item of trial t in run 1
        matching_trial = (tdata.imageNumber==data_R_s(:,2))&(subject_ID(s)==data_R_s(:,1));
        data_R_s(matching_trial,[4,7,12]) = [t-ntrials_R,tdata.response,tdata.onsetResponse-tdata.onsetDisplay];
    end
    for t = 1:ntrials_R %loop over trials of 3rd run
        tdata = cdata.ImageRating2(t);
        %find trial in run 3 that matches item of trial t in run 1
        matching_trial = tdata.imageNumber==data_R_s(:,2);
        data_R_s(matching_trial,[5,8,13]) = [t,tdata.response,tdata.onsetResponse-tdata.onsetDisplay];
    end
    data_R_s(:,9) = mean(data_R_s(:,6:7),2); %avg ratings for runs 1 and 2
    data_R_s(:,10) = mean(data_R_s(:,6:8),2); %avg ratings for all runs
    data_R((s-1)*ntrials_R+1:s*ntrials_R,:) = data_R_s; %fill in data for all subjects
    
	%%% --- %%% --- %%% --- %%%
    %%% --- Encoding (3 blocks * 4 repetitions * 16 items = 192 trials)
    %%% --- %%% --- %%% --- %%%
    data_E_s = zeros(ntrials_E,size(data_E,2)); %only for current subject
    image_map = cdata.Experiment.imageMap; %this is the link between image "identities" and their locations
    for r = 1:4 %loop over encoding runs
        for t = 1:ntrials_E
            T = (r-1)*ntrials_E+t; %"real" trial number (taking run into account)
            tdata = cdata.Encoding(T);
            image_blocknumber = tdata.gridBlockNumber;
            image_gridnumber = tdata.imageGridNumber;
            image_id = image_map((image_blocknumber-1)*nitems_per_block+image_gridnumber);
            if r == 1 %things are a bit different for the first run compared to the remaining runs
                data_E_s(t,1:5) = [subject_ID(s),image_id,image_blocknumber,image_gridnumber,t];
            else
                %find trial in current run that matches item of run 1
                matching_trial = image_id==data_E_s(:,2);
                data_E_s(matching_trial,3+r*2) = t;
            end
        end
        %proceed with encoding recall
        for t = 1:ntrials_E
            T = (r-1)*ntrials_E+t; %"real" trial number (taking run into account)
            tdata = cdata.EncodingRecall(T);
            image_blocknumber = tdata.gridBlockNumber;
            image_gridnumber = tdata.imageGridNumber;
            image_id = image_map((image_blocknumber-1)*nitems_per_block+image_gridnumber);
            rt_blank = tdata.onsetBlankResponse-tdata.onsetBlankDisplay;
            rt_image = tdata.onsetImageResponse-tdata.onsetImageDisplay;
            %find trial in current run that matches item of encoding run 1
            matching_trial = image_id==data_E_s(:,2);
            data_E_s(matching_trial,[4,11,12]+r*2) = [t,rt_blank,rt_image];
        end
    end
    data_E((s-1)*ntrials_E+1:s*ntrials_E,:) = data_E_s; %fill in data for all subjects
    
	%%% --- %%% --- %%% --- %%%
    %%% --- Decision (3 blocks * 40 decisions = 120 trials)
    %%% --- %%% --- %%% --- %%%
    data_D_s = zeros(ntrials_D,size(data_D,2)); %only for current subject
    control_block = cdata.Experiment.controlBlockNumber;
    for t = 1:ntrials_D %loop over trials
        tdata = cdata.Decision(t);
        %get which items have been presented
        image_blocknumber = tdata.gridBlockNumber;
        is_control_block = image_blocknumber==control_block; %whether items are shown in current trial 
        image_gridnumber_l = tdata.leftImageGridNumber;
        image_gridnumber_r = tdata.rightImageGridNumber;
        image_id_l = image_map((image_blocknumber-1)*nitems_per_block+image_gridnumber_l);
        image_id_r = image_map((image_blocknumber-1)*nitems_per_block+image_gridnumber_r);
        %get average rating values
        avg_ratings_l = data_R_s(data_R_s(:,2)==image_id_l,9:10);
        avg_ratings_r = data_R_s(data_R_s(:,2)==image_id_r,9:10);
        %get decision (recoded) and response time
        choice = find(ismember({'left','right','invalid'},tdata.choiceResponse));
        rt = tdata.onsetChoiceResponse-tdata.onsetChoiceDisplay;
        %put all together
        data_D_s(t,:) = [subject_ID(s),image_blocknumber,is_control_block,t,...
                         image_id_l,image_id_r,avg_ratings_l,avg_ratings_r,choice,rt];
    end
    data_D((s-1)*ntrials_D+1:s*ntrials_D,:) = data_D_s; %fill in data for all subjects

	%%% --- %%% --- %%% --- %%%
    %%% --- Cued Recall (3 blocks * 16 items = 48 trials)
    %%% --- %%% --- %%% --- %%%
    data_C_s = zeros(ntrials_C,size(data_C,2)); %only for current subject
    for t = 1:ntrials_C %loop over trials
        tdata = cdata.CuedRecall(t);
        %get which items have been presented
        image_blocknumber = tdata.gridBlockNumber;
        is_control_block = image_blocknumber==control_block; %whether items are shown in current trial 
        image_gridnumber = tdata.imageGridNumber;
        image_id = image_map((image_blocknumber-1)*nitems_per_block+image_gridnumber);
        %get behavioral information
        rt_blank = tdata.onsetBlankResponse-tdata.onsetBlankDisplay;
        confidence_rating = tdata.ratingResponse;
        rt_confidence = tdata.onsetRatingResponse-tdata.onsetRatingDisplay;
        %put all together
        data_C_s(t,:) = [subject_ID(s),image_blocknumber,is_control_block,t,...
                         image_id,rt_blank,confidence_rating,rt_confidence];
    end
    data_C((s-1)*ntrials_C+1:s*ntrials_C,:) = data_C_s; %fill in data for all subjects
    
    disp(['Done with subject ',int2str(s),' of ',int2str(nsubj)])
end

%%% --- %%% --- %%% --- %%%
%%% --- Generate csv files for R
%%% --- %%% --- %%% --- %%%

%rating
rating_titles = {'SubjectID','ImageNumber','TrialNumberR1','TrialNumberR2','TrialNumberR3','RatingR1','RatingR2',...
                 'RatingR3','AvgRatingR1R2','AvgRatingR1R2R3','RespTimeR1','RespTimeR2','RespTimeR3'};
table_R = array2table(data_R);
table_R.Properties.VariableNames(1:size(data_R,2)) = rating_titles;
writetable(table_R,'RATING.csv');

%encoding
encoding_titles = {'SubjectID','ImageNumber','BlockNumber','GridNumber','TrialNumberEnc1','TrialNumberEncRec1',...
                   'TrialNumberEnc2','TrialNumberEncRec2','TrialNumberEnc3','TrialNumberEncRec3','TrialNumberEnc4',...
                   'TrialNumberEncRec4','RespTimeBlank1','RespTimeImg1','RespTimeBlank2','RespTimeImg2',...
                   'RespTimeBlank3','RespTimeImg3','RespTimeBlank4','RespTimeImg4'};
table_E = array2table(data_E);
table_E.Properties.VariableNames(1:size(data_E,2)) = encoding_titles;
writetable(table_E,'ENCODING.csv');

%decision
decision_titles = {'SubjectID','BlockNumber','IsControlBlock','TrialNumber','LeftImageNumber','RightImageNumber',...
                   'LeftAvgRatingR1R2','LeftAvgRatingR1R2R3','RightAvgRatingR1R2','RightAvgRatingR1R2R3','Choice','RespTime'};
table_D = array2table(data_D);
table_D.Properties.VariableNames(1:size(data_D,2)) = decision_titles;
writetable(table_D,'DECISION.csv');

%cued recall
cued_recall_titles = {'SubjectID','BlockNumber','IsControlBlock','TrialNumber',...
                      'ImageNumber','RespTimeBlank','RecallRating','RespTimeRating'};
table_C = array2table(data_C);
table_C.Properties.VariableNames(1:size(data_C,2)) = cued_recall_titles;
writetable(table_C,'CUED_RECALL.csv');
    
end
