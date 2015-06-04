function progress(varargin)
% before loop t0=clock;
% in loop progress(sim_index,nSimulations,t0);

if nargin==2
    iteration=varargin{1};
    nIterations=varargin{2};
elseif nargin==3
    iteration=varargin{1};
    nIterations=varargin{2};
    startTime=varargin{3};
elseif nargin==4
    iteration=varargin{1};
    nIterations=varargin{2};
    startTime=varargin{3};
    userData=varargin{4};
else
    disp('incorrect number of inputs...')
end

if nIterations<8
    blockSize=round(nIterations/3);
else
    blockSize=round(nIterations/6);
end
if nargin>=2
    if mod(iteration,blockSize)==1
        fprintf('Progress: %04d / %04d (%03.2f%%)\n',[iteration nIterations iteration/nIterations*100])
        drawnow
    end
end

if nargin>=3
    if mod(iteration,blockSize)==1
        elapsed=etime(clock,startTime);
        perIteration=elapsed/iteration;
        remaining=perIteration*(nIterations-iteration);
        
        disp(['Time elapsed: ' sec2time(elapsed) ' - remaining: ' sec2time(floor(remaining)) ' (total: ' sec2time(round(elapsed+remaining)) ')'])
        drawnow
    end
end

if nargin>=4
    if mod(iteration,blockSize)==1
        disp('Extra info:')
        disp(sci(userData,2))
        drawnow
    end
end

