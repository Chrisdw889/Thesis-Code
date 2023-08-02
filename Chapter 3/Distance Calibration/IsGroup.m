%%% Version - v1.06
%%% Author  - Christopher Weaver

%%% Summary:
%%% IsGroup is a group study subclass that governs the processing of
%%% current-distance spectroscopy experiments. This includes I(t)
%%% experiments that utilise the spectroscopy window.

classdef IsGroup < GroupStudy
    
    properties
        time
        current %Column Vector
        distance
        traceNum
        tracesPoints
        servoSetPoint
    end
    
    methods
        
        function obj = IsGroup(directory, type)
            obj = obj@GroupStudy(directory);
            
            obj = obj.import(type);
        end
        
        
        function plotAllCurrentDistance(obj, figureNum)
            
            figure(figureNum);
            plot(obj.distance, obj.current);
        end
        
        function plotAllCurrentTime(obj, figureNum)
            
            figure(figureNum);
            plot(obj.time, obj.current);
        end
        
        
        function events = eventSearch(obj, baseRangeDev, eventTriggerDev, eventDurCutoff)
            
            startPoints = [];
            endPoints = [];
            traceIndices = [];
            
            for i = 1:obj.traceNum
                
                trace = obj.current(:, i);
                dev = std(trace);
                baseline = mean(trace) + (dev * baseRangeDev);
                eventTrigger = baseline + (dev * eventTriggerDev);
                
                j = 1;
                while j <= obj.tracesPoints
                    val = trace(j);
                    
                    if val > eventTrigger
                        
                        points = IsGroup.scopeEvent(trace, j, baseline, eventTrigger, obj.tracesPoints);
                        
                        
                        startPoint = points(1);
                        endPoint = points(2);
                        
                        
                        eventDurPoints = endPoint - startPoint;
                        if eventDurPoints < eventDurCutoff
                            continue;
                        end
                        
                        startPoints = [startPoints, startPoint];
                        endPoints = [endPoints, endPoint];
                        traceIndices = [traceIndices, i];
                        
                        j = endPoint + 1;
                        
                    else
                        j = j + 1;
                    end
                end
            end
            
            events = [startPoints', endPoints', traceIndices'];
            
        end
        
        
        
        
        function ret = betaAnalysis(obj)
            
            logCurrent = log(obj.current);
            figure();
            plot(obj.distance, logCurrent);
            
            %             thresholdParam = 0 % bool 0 = dist, 1 = current
            
            thresholdWord = input('Please input the desired threshold parameter. (Current or Distance)\n');
            
            NoiseThreshhold = 0.0;
            limit = 2;
            
            switch thresholdWord
                case 'Current'
                    NoiseThreshhold = input('Please input the current at which noise is measured. (in nA)\n');
                    NoiseThreshhold = NoiseThreshhold;
                    
                    val = obj.current(1, 1);
                    
                    while val > NoiseThreshhold
                        if limit == obj.tracesPoints
                            break;
                        end
                        %                         disp('limit = ' + limit);
                        val = obj.current(limit, 1);
                        limit = limit + 1;
                    end
                    
                case 'Distance'
                    NoiseThreshhold = input('Please input the distance at which noise is measured. (in nm)\n');
                    NoiseThreshhold = NoiseThreshhold * 1E-9;
                    val = obj.distance(1, 1);
                    
                    while val < NoiseThreshhold
                        if limit == obj.tracesPoints
                            break;
                        end
                        val = obj.distance(limit, 1);
                        limit = limit + 1;
                    end
                    
                    
            end
            
            figure();
            plot(obj.distance(1:limit, :), logCurrent(1:limit, :));
            
            
            
            
            betas = zeros(1, obj.traceNum);
            
            for i=1:obj.traceNum
                logCurr = logCurrent(1:limit, i);
                dist = obj.distance(1:limit, i) / 1E-9; % in nm
                
                tmp = logCurr;
                tmp = real(tmp);
                tmp(tmp == inf) = 2;
                tmp(tmp == -inf) = -2;
                logCurr = tmp;
                
                disp(i);
                f = fit(dist, logCurr, 'poly1');
                
                betas(1, i) = -f.p1;
            end
            
            figure();
            m_hist = histogram(betas, 10);
            
            figure();
            plot(betas);
            
            ret = betas;
            
        end
        
        
        % 0 Break Junctions
        % 1 Time Spec
        function obj = import(obj, fileType)
            
            header = "";
            
            if fileType == 0
                header = "I(s)*.txt";
            elseif fileType == 1
                header = "I(t) Spec*.txt";
            end
            
            home = "Y:\CDW\Thesis Work by Chapter\I(t) Classification\Distance Calib\MATLAB Scripts";
            
            cd(obj.directory);
            files = dir(header);
            cd(home);
            
            filesNum = length(files);
            k = 0; % TraceNumber
            
            pointsPerTrace = 0;
            
            for i=1:filesNum
                tmpData = obj.importFile(files(i).name);
                dataStartRow = 0;
                chunks = 0;
                
                for j=1:length(tmpData(:, 1))
                    entry = tmpData(j, 1);
                    
                    if contains(entry, "servoSetpoint")
                        spEntry = tmpData(j, 2);
                        %                         splitEntry = split(entry, " ");
                        %                         len2 = length(splitEntry);
                        obj.servoSetPoint = str2double(spEntry);
                    end
                    if entry == "chunk"
                        chunks = chunks + 1;
                        pointsPerTrace = str2double(tmpData(j, 3));
                    end
                    if entry == "data"
                        dataStartRow = j + 2;
                        break;
                    end
                    
                end
                
                allTime = str2double(tmpData(dataStartRow:end, 1));
                allDistance = str2double(tmpData(dataStartRow:end, 2));
                allCurrent = str2double(tmpData(dataStartRow:end, 3));
                
                timeMat = reshape(allTime, [pointsPerTrace, chunks]);
                distanceMat = reshape(allDistance, [pointsPerTrace, chunks]);
                currentMat = reshape(allCurrent, [pointsPerTrace, chunks]);
                
                obj.time = horzcat(obj.time, timeMat);
                obj.distance = horzcat(obj.distance, distanceMat);
                obj.current = horzcat(obj.current, currentMat);
                
                k = k + chunks;
                
            end
            
            obj.tracesPoints = pointsPerTrace;
            obj.traceNum = k;
            
            
        end
        
        function data = importFile(obj, filename, startRow, endRow)
            delimiter = {'\t',' '};
            
            if nargin<=2
                startRow = 1;
                endRow = inf;
            end
            
            formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';
            
            fileID = fopen(obj.directory + "\" + filename,'r');
            
            dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            for block=2:length(startRow)
                frewind(fileID);
                dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
                for col=1:length(dataArray)
                    dataArray{col} = [dataArray{col};dataArrayBlock{col}];
                end
            end
            
            fclose(fileID);
            
            data = [dataArray{1:end-1}];
        end
    end
    
    methods(Static)
        
        function points = scopeEvent(currenttrace, k, baseline, eventTrigger, tracepoints)
            
            l = k;
            m = k;
            while l > 1
                if currenttrace(l) < baseline
                    % Search back 5
                    count = 5;
                    restart = false;
                    index = l - 1;
                    while index > 0 && count >= 0
                        val2 = currenttrace(index);
                        if val2 > eventTrigger
                            restart = true;
                            break;
                        end
                        index = index - 1;
                        count = count - 1;
                    end
                    if restart
                        l = index;
                    else
                        break;
                    end
                else
                    l = l - 1;
                end
                
                
            end
            
            while m < tracepoints
                if currenttrace(m) < baseline
                    count = 5;
                    restart = false;
                    index = m + 1;
                    while index <= tracepoints && count >= 0
                        val2 = currenttrace(index);
                        if val2 > eventTrigger
                            restart = true;
                            break;
                        end
                        index = index + 1;
                        count = count - 1;
                    end
                    if restart
                        m = index;
                    else
                        break;
                    end
                else
                    m = m + 1;
                end
            end
            
            points = [l, m];
            
        end
        
    end
    
    
end