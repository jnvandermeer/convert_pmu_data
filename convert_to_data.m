

% pre-requisites
% eeglab - WITH fmrib plugin
% spm
% bash
% a working matlab on Linux platform


TMPFILENAME='tmp.txt';

% ask user for the physiological file
% ask where the ap fieldmap is
[pl ~] = spm_select(1,'any','Select a physiological logfile',{},pwd,'.*');
[plp plb ple] = fileparts(pl);

dfiles=dir([plp filesep plb '.*']);
fprintf('found %d files \n',numel(dfiles));

[dd ~] = spm_select(1,'dir','Select directory of associated DICOM IMAs',{},pwd,'.*');

fprintf('running bash command to obtain some timestamp info\n');
s=['for f in `ls ' dd filesep '*.IMA`;do echo dcminfo $f -tag 0008 0032;done |sh > ' TMPFILENAME];
system(s);
fid=fopen(TMPFILENAME);
lines={};
while ~feof(fid)
    lines{end+1}=fgetl(fid);
end
fclose(fid);
delete(TMPFILENAME);

ts_dicoms=[];
for i=1:numel(lines)
    ts_dicoms(end+1) =  str2num(lines{i}(1:2)) * 3600 * 1000 + ...
                        str2num(lines{i}(3:4)) * 60 * 1000 + ...
                        str2num(lines{i}(5:6)) * 1000 + ...
                        floor(str2num(lines{i}(8:10)));
end




datavecs={};
datanames={};

for i=1:numel(dfiles);

    [~, fn,ext] = fileparts(dfiles(i).name);

    fname = [plp  filesep fn ext];
    
    if ~strcmp(ext,'.pmu')&&~strcmp(ext,'.set');
        
        
        
        fid=fopen(fname);
        
        lines={};
        while ~feof(fid);lines{end+1}=fgetl(fid);end;fclose(fid);

        % miliseconds since the New Day:
        ts_start=str2num(regexprep(lines{13},'LogStartMDHTime:\s*',''));
        ts_stop =str2num(regexprep(lines{14},'LogStopMDHTime:\s*',''));
        
        fprintf('timestamp info...: %d %d \n',ts_start,ts_stop)
        
        
        clear s;
        s=regexp(lines{1},'[^\s]*','match');

        
        d=[];
        fprintf('reading data from file %s\n',fname);
        
        % because it starts with 7th number!!
        for j=7:numel(s);d(end+1)=str2double(s{j});end
        
        
        if strcmp(ext,'.ecg')
        % keyboard;
        end
        % getting rid of triggers...
        excludes=[find(d==3785) find(d==5000) find(d==5003) find(d==6000) find(d==6002)];
        d(excludes)=[];
        
        
        
        if strcmp(ext,'.ecg')
            
            
            r=rem(numel(d),4);
            if r>0
                fprintf('getting rid of %d final samples to fit into 4 rows\n',r);
                toberemoved=(numel(d)-r+1):numel(d);
            else
                toberemoved=[];
            end    
            d(toberemoved)=[];
            
            
            m=reshape(d,4,numel(d)/4)';
            
            % apply very dirty trick
            % we will ... SORT!!!
            % my first steps into the dark side is complete.
            % and i also use ugly k
            
            
            % THEN - sort!
            m2=[];
            for k=1:size(m,1)
                m2(k,:) = sort(m(k,:));
            end
            m=m2;
            
            for k=1:size(m,2)
                datavecs{end+1}=m(:,k)';
                datanames{end+1}=sprintf('ecg%d',k);
            end
            
            
            % keyboard;
            
        else
            
            m=d;
            
            datavecs{end+1}=m;
            datanames{end+1}=ext(2:end);
            
        end
        
        
     
    else
        
        fprintf('.pmu or .set file! -- skipping...');
    end
    
    
end

% logfiles start at different times but stop at at least hte same time.

% Sooo.. let's figure out our sync while counting backwards from ts_stop.
% taking into account 400 Hz sampling rate.


for i=1:numel(datavecs)
    datavecs{i}(isnan(datavecs{i}))=mean(datavecs{i});
end

% figure out how much data we can manage - while putting it into one huge
% array.
% ALSO... the ECG starts at a different time from the resp at a different
% time from the PPU. So.. this also makes sense, even.
t=[];for i=1:numel(datavecs);t(end+1)=numel(datavecs{i});end
datasize=min(t);

% initialize eeglab..
% make an empty dataset...
EEG=eeg_emptyset();

for i=1:numel(datanames);
    EEG.chanlocs(i).labels=datanames{i};
end
    
EEG.data=[];
for i=1:numel(datavecs)
    EEG.data(i,:)=datavecs{i}(end-datasize+1:end);
end
EEG.pnts=size(EEG.data,2);
EEG.srate=400;
EEG.setname=plb;

% keyboard;


% NOW -- markers of dicoms.
EEG.event=struct();
EEG.event(:)=[];
for i=1:numel(ts_dicoms)
    
    EEG.event(end+1).type='mri';
    EEG.event(end).latency=EEG.pnts-floor((ts_stop-ts_dicoms(i)+mean(diff(ts_dicoms)))/1000*400);
end

EEG=eeg_checkset(EEG);
%EEG = pop_eegfiltnew(EEG, [],0.05,26400,1,[],1);

for i=1:EEG.nbchan
    v=EEG.data(i,:);
    v2=v;v2(isnan(v2))=[];
    v(isnan(v))=mean(v2);
    EEG.data(i,:)=v;
end
    
EEG.data=detrend(EEG.data','linear')';
EEG=eeg_checkset(EEG);


% NOW... correct the ECG:
EEG2 = pop_select( EEG,'channel',{'ecg1' 'ecg2' 'ecg3' 'ecg4'});
EEG2 = eeg_checkset( EEG2 );
EEG2 = pop_fmrib_fastr(EEG2,70,10,30,'mri',0,1,0,0,0,0.03,[],'auto');
EEG2 = eeg_checkset( EEG2 );

% and replace the ECG:
EEG.data(1:4,:)=EEG2.data(1:4,:);


% and do some detection(s) on the ECG channel that's best (by default,
% that's # 4)
EEG=pop_fmrib_qrsdetect(EEG,4,'qrs','no');
EEG = eeg_checkset( EEG );
pop_eegplot( EEG, 1, 1, 1);



EEG=pop_saveset(EEG,[plp filesep plb '.set']);


