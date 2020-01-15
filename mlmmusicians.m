%{
    file_name : mlmmusicians.m
    author : Francesco Mantegna
    institution : NYU
    project : Music&Poetry
    date : 12/01/2019
%}

%% add fieldtrip path

full_path_to_fieldtrip = '/Users/francesco/Documents/MATLAB/fieldtrip-20170414';
addpath(full_path_to_fieldtrip)
ft_defaults

%% defining variables

twin      =    [-0.50 1.0];
rawDir    =    '/Users/francesco/Documents/MATLAB/MPI/MusicPrediction/RawMusicians';
inputDir  =    '/Users/francesco/Documents/MATLAB/MPI/MusicPrediction/DataMusic/Data0.1-30MUSICIANS';
excelDir  =    '/Users/francesco/Documents/MATLAB/MPI/MusicPrediction/mlm_music';
numlist   =    [1, 3, 4, 8, 11, 12, 14, 15, 16, 17, 18, 19, 21, 32, 34, 35, 36, 37]; 
subjnum   =    18; 

for s = 1:subjnum
    
    cd(rawDir)

    %% loading raw files

    cfg = [];
    if s <= 4 
        cfg.headerfile   = ['musicpred_m0' int2str(numlist(s)) '_s2.vhdr']; 
        cfg.datafile     = ['musicpred_m0' int2str(numlist(s)) '_s2.eeg'];
    elseif s >= 5 && s <= 15
        cfg.headerfile   = ['musicpred_m' int2str(numlist(s)) '_s2.vhdr'];
        cfg.datafile     = ['musicpred_m' int2str(numlist(s)) '_s2.eeg'];
    elseif s >= 16
        cfg.headerfile   = ['musicpred_m' int2str(numlist(s)) '.vhdr']; 
        cfg.datafile     = ['musicpred_m' int2str(numlist(s)) '.eeg'];
    end

    %% defining trials

    cfg.trialfun     = 'francesco_trialfun_ERPmusic';

    cfg.trialdef.pre  = abs(twin(1));
    cfg.trialdef.post = twin(2);
    cfg = ft_definetrial(cfg);
    if s == 11
        cfg.trl(41,:) = [];
        cfg.trl(47,:) = [];
    end
    
    rawdata = ft_preprocessing(cfg);
    
    % load pre-processed data
    
    cd(inputDir)
    
    load(['prepromusic_SUBJ', num2str(s),'.mat'])
    
    for pp = 1:length(data_no_artifacts.trial)
        pt = rawdata.sampleinfo(:,1) == data_no_artifacts.sampleinfo(pp,1);
        idx = find(pt);
        data_no_artifacts.trialinfo(pp,2) = rawdata.cfg.trl(idx,5) - 30;
    end
    
    data = data_no_artifacts;
    
    clear data_no_artifacts rawdata
    
    load('standard_1020_elec.mat')
    r = load('musiciansrating.mat');
    avgrating = r.musicians;
    x  = num2cell(elec.pnt(:,1));
    y  = num2cell(elec.pnt(:,2));
    z  = num2cell(elec.pnt(:,3));
    ch = char(elec.label);
    cellch = cellstr(ch);
    
    subj_id = num2cell(s + 31);
    group = {char('musicians')};
    
    for t = 1:length(data.trial)
        temp = cell(length(ch),12); 
        for c = 1:length(ch)
            n5 = num2cell(mean(data.trial{t}(c,476:551),2)); % N500 window 450-600 ms
            p3 = num2cell(mean(data.trial{t}(c,351:401),2)); % P300 window 200-300 ms
            bs = num2cell(mean(data.trial{t}(c,151:251),2)); % baseline window -250-0 ms
            if data.trialinfo(t,1) == 1
                condition = {char('tonic')};
                rating    = num2cell(avgrating(s,1));
            elseif data.trialinfo(t,1) == 2
                condition = {char('dominant')};
                rating    = num2cell(avgrating(s,2));
            else
                condition = {char('aug4')};
                rating    = num2cell(avgrating(s,3));
            end
            scale = num2cell(data.trialinfo(t,2));
              temp(c,:) = [scale, subj_id, group, condition, n5, p3, bs, x(c), y(c), z(c), cellch(c), rating];
        end
        if t == 1
            info = temp;
        else
            info = vertcat(info, temp);
        end
    end
    
    sortedinfo =  sortrows(info,1);
    if s == 1
        inputm = info;
    else
        inputm = vertcat(inputm, info);
    end

end

cd(excelDir)

save('mlminputm.mat', 'inputm')