classdef Import
    
    methods(Static)
        
        %% Functions
        function strucdata = import(directory, filename, filetype)
            
            datatype = -1; % set to 1 BINARY or 2 (ascii)
            
            fmt = Import.scopefile(directory, filename);
            datatype = fmt.datatype;
            dataline = fmt.datastart;
            
            header = Import.importheader(directory, filename, filetype, dataline);
            cols = length(header.columns);
            rows = header.NumPoints;
            
            
            data = [];
            if datatype == 1
                data = Import.importbin(directory, filename, filetype, dataline, cols, rows);
            elseif datatype == 2
                data = importasc();
            else
                disp("error");
            end
            
            %Rearrange data based on chunks
            numTraces = length(header.chunks);
            pointsPerTrace = header.chunks{1};
            
            shapeddata = reshape(data, pointsPerTrace, numTraces);
            
            strucdata.header = header;
            strucdata.data = shapeddata;
            strucdata.numTraces = numTraces;
            strucdata.pointsPerTrace = pointsPerTrace;
            
            
        end
        
        
        %% Helper Functions
        
        
        function fmt = scopefile(directory, filename)
            
            fid = fopen(strcat(directory, filename));
            
            tline = fgetl(fid);
            i = 1;
            while ischar(tline)
                tline = fgetl(fid);
                i = i + 1;
                
                if contains(tline, 'BINARY')
                    disp(tline)
                    datatype = 1;
                    break;
                    
                end
                if contains(tline, 'ASCII')
                    disp(tline)
                    datatype = 2;
                    break;
                end
                
            end
            
            fclose(fid);
            
            fmt.datastart = i + 1;
            fmt.datatype = datatype;
            
            
        end
        
        function header = importheader(directory, filename, filetype, dataline)
            fid = fopen(strcat(directory, filename));
            
            header.columns = cell(1);
            
            chunks = {};
            
            tline = fgets(fid);
            i = 1;
            while i < dataline
                tline = fgetl(fid);
                i = i + 1;
                
                if strcmp(filetype, 'CV')
                    
                elseif strcmp(filetype, 'I(t) Spec')
                    if contains(tline, "servoSetpoint")
                        header.setpoint = Import.extractHeaderVal(tline, "servoSetpoint");
                    end
                elseif strcmp(filetype, 'I(s)')
                    if contains(tline, "bias")
                        header.bias = Import.extractHeaderVal(tline, "bias");
                    end
                    if contains(tline, "servoSetpoint")
                        header.setpoint = Import.extractHeaderVal(tline, "servoSetpoint");
                    end
                else
                    disp("Unknown filetype");
                end
                
                if contains(tline, 'bufferLabel')
                    j = length(header.columns) + 1;
                    header.columns{j} = strrep(tline, 'bufferLabel', '');
                end
                
                
                if contains(tline, 'chunk')
                    tline = strrep(tline, ' ', '');
                    tline = strrep(tline, 'chunk', '');
                    
                    sections = strsplit(tline, '\t');
                    chunkindex = str2num(sections{1});
                    chunks{chunkindex + 1} = str2num(sections{2});
                    
                end
                
            end
            
            fclose(fid);
            
            header.NumPoints = sum(cell2mat(chunks));
            header.chunks = chunks;
        end
        
        function datamat = importbin(directory, filename, filetype, dataline, cols, rows)
            fid = fopen(strcat(directory, filename));
            
            %Column correction
            if strcmp(filetype, 'I(s)')
                cols = cols - 2;
            elseif strcmp(filetype, 'I(t) Spec')
                cols = cols - 2;
            end
            
            i = 1;
            while i < dataline
                tline = fgetl(fid);
                i = i + 1;
            end
            
            datamat = fread(fid, 'float32');
            
            fclose(fid);
            
            
        end
        
        function val = extractHeaderVal(txtline, key)
%             tline = strrep(txtline, ' ', '\t');
%             tline = strrep(tline, key, '');
%             valstr = strsplit(tline, '\t');
            val = sscanf(txtline, strcat(key,'%e'));
        end
    end
    
end