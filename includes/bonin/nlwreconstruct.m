function rec = nlwreconstruct(stack,shift, dsratio)
%NLWRECONSTRUCT
% REC = NLWRECONSTRUCT(STACK)
% REC = NLWRECONSTRUCT(STACK,SHIFT)
% REC = NLWRECONSTRUCT(STACK,SHIFT, DSRATIO)

if nargin < 2
    shift = 0;
end;

if nargin < 3 
    dsratio = 1;
end;

if isempty(intersect(dsratio,[1 2 4]))
    error('dsratio should be 1, 2 or 4');
end;

[nCycles,nSamples,nFrames]=size(stack);

%%
switch nSamples % # samples per cyle
    case 9000  % 80 MHz   
       lineWidth = 2880/dsratio; % chosen so no sample discarded at FOV center
    case 2308  % 20 MHz
        lineWidth = 720/dsratio; % chosen so no sample discarded at FOV center
    otherwise
        nSamples % added output for debug
        error('unknown raw stack dimensions');
end

period = nSamples; % # samples per cycle
delta = 2*pi / period; % phase increment per sample
phi = [0:round(period/2)-1]*delta; % scanner phase


%%
rec = int16(zeros(nCycles*2,lineWidth,nFrames));

% format input parameters
binning_tbl = uint16( (1 - cos(phi))*(lineWidth-1)/2 + 1 ); % range 1:lineWidth
% plot(diff(binning_tbl),'.-k')
line_start_ind = uint32([0:nCycles-1]*nSamples+1);
nSamplesPerLine = uint32(period/2);
lineWidthU = uint16(lineWidth);
delayCorrect = int32(shift);

tic;
for iFr = 1:nFrames
    frame = stack(:,:,iFr)';
    data = (frame(:));
    
    %rec(:,:,iFr) = mymex(data,...
    rec(:,:,iFr) = nlwformatframe(data,...
                               line_start_ind,...
                               binning_tbl,...
                               lineWidthU,...
                               nSamplesPerLine,...
                               delayCorrect);
end;
t = toc;
fprintf(1,'%2.1f fps \n',nFrames/t)
% imagesc(mean(rec,3));

return;
