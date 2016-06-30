classdef log_object < handle
    properties
        entry=struct('name','','time',[],'elapsed',[],'options',[]);
        nEntries=0;        
    end
    
    methods
        function self=log_object(varargin)
            disp('Started new log file!')
        end
        
        
        function append(self,varargin)
            N=self.nEntries+1;
                        
            if nargin>=2
                self.entry(N).name=varargin{1};                
            end
            self.entry(N).time=clock;
            if nargin>=3
                self.entry(N).elapsed=etime(clock,varargin{2});
            end
            if nargin>=4
                self.entry(N).options=varargin{3};
            end
            self.nEntries=N;
        end
    end
end
