%{
    file_name : mlmcsv.m
    author : Francesco Mantegna
    institution : NYU
    project : Music&Poetry
    date : 12/01/2019
%}

excelDir  =    '/Users/francesco/Documents/MATLAB/MPI/MusicPrediction/mlm_music';
cd(excelDir)

load('mlminputm.mat')
load('mlminputnnm.mat')

titles = {'scale_id','subject_id','group','condition','N5', 'BS', 'x', 'y', 'z', 'channel', 'rating'};

mlminput = vertcat(inputnnm, inputm);
mlminput = sortrows(mlminput,2);
T = cell2table(mlminput);
T.Properties.VariableNames = titles;
writetable(T, 'mlm_inputmusic.csv');
