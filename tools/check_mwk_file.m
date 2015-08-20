clear all
clc

mwk_file_names={...
    '/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-08-14_AH03/2015-08-14_AH03.mwk',...
    '/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-08-07_AH03/2015-08-07_AH03.mwk',...
    '/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-08-19_AH02/2015-08-19_AH02.mwk',...
    '/Users/benvermaercke/Dropbox (coxlab)/2p-data/2015-08-19_AH03/2015-08-19_AH03.mwk'};

nFiles=length(mwk_file_names);

for iFile=1:nFiles
    mwk_file_name=mwk_file_names{iFile};
    A=getCodecs(mwk_file_name);
    event_codec=A.codec;
    %{event_codec.tagname}'
    
    E=get_events_by_name(mwk_file_name,'#stimDisplayUpdate',event_codec);
    %%
    T=double(cat(1,E.time_us));
    T_trial=diff(T)/1e6;
    T_trial(~between(T_trial,[0.5 1.5]))=[];
    %T_trial(~between(T_trial,[1.5 2.5]))=[];
    
    delta_frames=T_trial/(1/60);
    %delta_frames=delta_frames-min(delta_frames);
    subplot(nFiles,1,iFile)
    %plot(delta_frames)
    refresh_range=60:70;
    %refresh_range=120:130;
    hist(delta_frames,refresh_range)
    mwk_file_name_title=mwk_file_name;
    mwk_file_name_title=strrep(mwk_file_name_title,'_',' ');
    title(mwk_file_name_title)
    xlabel(sprintf('N=%d (std=%3.4f)',[length(T_trial) std(T_trial)]))
    box off
    axis([refresh_range([1 end]) 0 3000])
end
