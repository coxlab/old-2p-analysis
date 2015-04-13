classdef HandleObject < handle
   properties
      Object=[];
   end
 
   methods
      function obj=HandleObject(receivedObject)
         obj.Object=receivedObject;
      end
   end
end

% source: http://www.matlabtips.com/how-to-point-at-in-matlab/