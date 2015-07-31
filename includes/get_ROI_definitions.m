function varargout=get_ROI_definitions(varargin)

self=varargin{1};

%%% Newer data format allows for multiple ROI_definitions, second input
%%% argument can define which is selected. Default 1, value is set in
%%% header_script
if nargin>=2&&~isempty(varargin{2})   
    ROI_definition_nr=varargin{2};
else
    ROI_definition_nr=self.ROI_definition_nr;
end

if isfield(self,'ROI_definitions')
    if isfield(self.ROI_definitions,'ROI_nr')
        %disp('Reading old data format')
        ROIs=self.ROI_definitions;
    end
    
    if isfield(self.ROI_definitions,'ROI')
        %disp('Reading new data format')
        ROIs=self.ROI_definitions(ROI_definition_nr).ROI;
    end
    
    varargout{1}=ROIs;
elseif isprop(self,'ROI_definitions')
    if isfield(self.ROI_definitions,'ROI')
        ROIs=self.ROI_definitions(ROI_definition_nr).ROI;
        
        if ~isempty(ROIs(1).ROI_nr)
            varargout{1}=ROIs;
        else
            varargout{1}=[];
        end
    else
        varargout{1}=[];
    end    
else
    error('No ROIs field found...')
end