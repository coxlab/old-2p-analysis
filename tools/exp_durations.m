clear all
clc

expType=2;

switch expType
    case 1 % Leuven shapes
        stim_durations=[1 2]; % seconds, stim and blank
        nRows=3;
        nCols=3;
        nStimuli=6;
        nRepeats=10;
    case 2 % RSVP
        stim_durations=[1 2]; % seconds, stim and blank
        nRows=4;
        nCols=8;
        nStimuli=12;
        nRepeats=3;
end
nPositions=numel(meshgrid(1:nRows,1:nCols));

exp_duration_seconds=sum(stim_durations)*nPositions*nStimuli*nRepeats;
exp_duration_minutes=exp_duration_seconds/60;

nFrames=1500; % scanimage frames
frame_rate=3; % Hz
block_duration=nFrames/frame_rate;
nBlocks=ceil(exp_duration_seconds/block_duration);

fprintf('Presenting %d repeats for %d stimuli \nlasting %dseconds (%dsON/%dsOFF) \nat %d positions,\nwill take %3.2f minutes (%d sec).\nWhich will require at least %d blocks of %d frames at %dHz.\n',[nRepeats nStimuli sum(stim_durations) stim_durations nPositions exp_duration_minutes exp_duration_seconds nBlocks nFrames frame_rate])
