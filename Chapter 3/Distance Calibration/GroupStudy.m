%%% Version - v1.00
%%% Author  - Christopher Weaver

%%% Summary:
%%% GroupStudy is an intermediate abstract class that is to be inherited by
%%% studies that analyse multiple experiment files per object

classdef GroupStudy < Study
    
    properties
        directory
    end
    
    methods
        
        function obj = GroupStudy(directory)
            obj.directory = directory;
        end
      
    end
    
    
end