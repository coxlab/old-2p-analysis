%%% this script will convert .mwk data files to .mat files, everything or parts of it
%%% only works on mac systems so far.
%%% perhaps, we could build a python version that is portable across
%%% systems if python can export to .mat files

%%% in principle, we could crawl the entire data directory and convert
%%% every .mwk file we encounter. perhaps this will be the default behavior
%%% if just folder en not file is provided as input.

%%% We are in the process of trying to compile the mwfeval function for
%%% windows and linux, no luck so far... code seems to be outdated

clear all
clc

if ismac
    overwrite_existing=0;
    data_folder='/Users/benvermaercke/Dropbox (coxlab)/2p-data';
    mwk_file_name_list=rdir([data_folder '/**/**.mwk']);
    nFiles=length(mwk_file_name_list);
    for iFile=1:nFiles
        mwk_file_name=mwk_file_name_list(iFile).name;
        
        disp(['Processing: ' mwk_file_name])
        
        ok=1;
        if length(strfind(mwk_file_name,'.mwk'))==2
            mwk_file_name=fileparts(mwk_file_name);
        elseif length(strfind(mwk_file_name,'.mwk'))>2
            % skip when mwk is already processed more than once
            ok=0;
        end
        
        mat_file_name=strrep(mwk_file_name,'mwk','mat');
        if overwrite_existing==0&&exist(mat_file_name,'file')
            % skip when mat file exists
            ok=0;
        end
        
        if ok==1                        
            codecs=getCodecs(mwk_file_name);
            
            tag_names={codecs.codec.tagname}';
            codes=[codecs.codec.code]';
            
            event_code_strings={'#stimDisplayUpdate','ExpType','ExpName_short','stm_pos_x','show_vertical_bar'};
            [a,b]=ismember(event_code_strings,tag_names);
            b(b==0)=[];
            event_code_selection=codes(b);
            
            events=getEvents(mwk_file_name,event_code_selection);
            
            save(mat_file_name,'codecs','event_code_strings','event_code_selection','events')
            disp(['Saved to : ' mat_file_name])
        else
            disp('Skipped...')
        end
    end
else
    disp('Currently, this script only run on mac systems given that the mex file needed to read the .mwk file is only compiled for mac')
end