%% Log of code

% 23 Feb 2021, modified to have output as the start and duration of each
% sleep staging;

% 23 Feb 2021, modified to extract all arousal events;


%%
classdef loadPSGAnnotationClass 
    
    %---------------------------------------------------- Public Properties
    properties (Access = public)
        % Input
        fileName = ''; % eg: '123.edf'   
        vendorName = ''; % eg: 'Embla'
        sleepStageValues = [];
        sleepSpindleValues = [];%yh add spindle info 10/16/17
        annotationType = ''; % eg: 'PSGAnnotation', 'CMPStudyConfig'
        isSDO = 0;
        
        % Optional Parameters
        errMsg = {};
        % Error list
        errList = {};
        errMap = containers.Map('KeyType', 'char', 'ValueType', 'char');        
    end
    %------------------------------------------------- Dependent Properties
    properties (Dependent = true)        
        % PhysioMiMi Terms
        xmlEntries                        % XMl entry types
        ScoredEvent                       % Scored event structure list
        EventList                         % List of events
        EventTypes                        % Unique event entry types
        EventStart                        % Event start list
        SleepStages
        EpochLength
        SleepStagesStart
        SleepStagesDuration
        %mod2
        ArousalsStart
        ArousalsDuration
    end
    %------------------------------------------------- Protected Properties
    properties (Access = protected)  
        AnnotationType = ''; % eg: 'PSGAnnotation'
        EventConcepts = [];
        EventStages = [];
        
        % Lights off/on text
        lightsOffText = 'Lights Off';
        lightsOnText = 'Lights On';
        
        % PhysioMiMi Terms
        xmlEntriesP                   % XML Doc Entries
        ScoredEventP                  % Scored event data structure list
        EventListP                    % List of events names(EventConcepts)
        EventTypesP                   % Unique event names(categories)
        EventCategoryListP            % Event categories
        EventCategoryP                % Unique event categories
        EventStartP                   % Event start list
        SleepStagesP                  % Sleep Stages list
        EpochLengthP                  % Epoch Length  
        SleepStagesSt
        SleepStagesDur
        %mod2
        ArousalSt
        ArousalDur
    end
    %------------------------------------------------------- Public Methods
    methods
        %------------------------------------------------------ Constructor
        function obj = loadPSGAnnotationClass(varargin)
            if nargin == 1
                obj.fileName = varargin{1};
            else
                fprintf('obj = loadPSGAnnotationClass (filename)') % (fileName) to (filename)
            end 
        end
        
        %---------------------------------------------------------
        %available event names for public usage
        function eventNameList = availableEventNames(obj)
            % remove stages:
            eventNameList = obj.EventTypesP; %%% Added unique func, TODO
        end
        function eventTypeList = availableEventTypes(obj)
            eventTypeList = obj.EventCategoryP;
        end
        %--------------------------------------------------------- loadFile
        function obj = loadFile(obj)
            % Load Grass file
            filename = obj.fileName;
            fid = fopen(filename);      
            % Process if file is open
            if fid > 0
                fileTxt = fread(fid)';
            else
                msg = sprintf('Could not open %s', filename);
                error(msg);
            end
            %-----------------------------------------------  Resolve 'SDO'
            % Temp = fread(fid,[1 inf],'uint8'); % not efficient
            obj.isSDO = strfind(fileTxt,'SDO:');
            %-----------------------------------------------------
            
            % Pass loaded information to object
            try
                xdoc = xmlread(filename);
            catch
                errMsg{end+1} = 'Failed to read XML file';
                error('Failed to read XML file %s.',xmlfile);
            end
            [ScoredEvent, SleepStageNames, EpochLength, obj.sleepStageValues, obj.sleepSpindleValues,SleepStageDur, SleepStageSt,ArousalSt, ArousalDur] = ...
                parseAndValidateNodes(xdoc);            
           
            % Get event information
            eventListF = @(x)ScoredEvent(x).EventConcept;
            EventListP = arrayfun(eventListF,[1:length(ScoredEvent)],...
                'UniformOutput', 0)';
            EventTypesP = unique(EventListP);
            
            eventCategoryF = @(x)ScoredEvent(x).EventType;
            EventCategoryListP = arrayfun(eventCategoryF,[1:length(ScoredEvent)],...
                'UniformOutput', 0)';
            EventCategoryP = unique(EventCategoryListP);
            
            eventStartF = @(x)ScoredEvent(x).Start;
            EventStartP = arrayfun(eventStartF,[1:length(ScoredEvent)],...
                'UniformOutput', 0)';

            % Pass key varaibles to obj
            obj.ScoredEventP = ScoredEvent;
            obj.SleepStagesP = SleepStageNames';
            obj.EpochLengthP  = EpochLength;
            obj.SleepStagesDur = SleepStageDur;
            obj.SleepStagesSt = SleepStageSt;
            % Mod2
            obj.ArousalDur = ArousalDur;
            obj.ArousalSt = ArousalSt;
            
            % Pass detail information to obj
            obj.EventListP = EventListP;
            obj.EventTypesP = EventTypesP;
            obj.EventStartP = EventStartP;
            obj.EventCategoryP = EventCategoryP;
            
            fid = fclose(fid);
        
            %------------------------------------ Parse and validate nodes
            function [ScoredEvent, SleepStageNames, EpochLength, SleepStageValues, SleepSpindleValues, SleepStageDur,SleepStageSt, ArousalSt, ArousalDur] = ...
                    parseAndValidateNodes(xmldoc)
                fprintf('\n>>> Parsing annotation file... \n')
                % Function parses each XML node
                xmlVersion = xmldoc.getXmlVersion;
                xmlEncoding = xmldoc.getXmlEncoding;
                rootNode = xmldoc.getFirstChild;
                rootNodeTag = rootNode.getTagName;
                obj.AnnotationType = rootNodeTag;

                Temp = xmldoc.getElementsByTagName('EpochLength');
                EpochLength = str2double(Temp.item(0).getTextContent);
                
                TempVendor = xmldoc.getElementsByTagName('SoftwareVersion');
                obj.vendorName = TempVendor.item(0).getTextContent;                               
                
                events = xmldoc.getElementsByTagName('ScoredEvent');
                % Add check code to deal with the missing subfield(Start/Duration, etc)
                if events.getLength > 0
                    SleepStageNames = {};
                    SleepStageDur = {};
                    ArousalSt = {};
                    ArousalDur = {};
                    SleepStageSt = {};
                    SleepStageValues = [];
                    SleepSpindleValues = [];
                    ScoredEvent = struct(...
                        'EventType', '', ...
                        'EventConcept', '', ...
                        'Start', [], ...
                        'Duration', [], ...
                        'InputCh', '', ...
                        'SpO2Baseline', [], ...
                        'SpO2Nadir', [], ...
                        'ToolTip', '' ...
                    );

                    stagesNameVector = {...
                        'Stage 1 sleep',...
                        'Stage 2 sleep',...
                        'Stage 3 sleep',...
                        'Stage 4 sleep',...
                        'REM sleep',...
                        'Wake',...
                        'Movement'
                    };
                    
                    for i = 0 : events.getLength - 1
                        
                        eventConceptText = '';
                        nadirNum = [];
                        baselineNum = [];
                        startNum = [];
                        durationNum = [];
                        hasDesaturation = [];
                        eventValid = 1;
                        InputChName = '';
                        %%% First check if the event has predefined event
                        %%% length, if not report error(May need a report event error mechanism)
                        try
                            eventConceptNode = events.item(i).getElementsByTagName('EventConcept');
                            eventConceptText = char(eventConceptNode.item(0).getTextContent);
                            tooltip = eventConceptText;
                            if i ~= 0
                                tmp = strsplit(eventConceptText, '|');
                                eventConceptText = tmp{2};
                                tooltip = tmp{1};
                            end
                            
                            eventTypeNode = events.item(i).getElementsByTagName('EventType');
                            eventTypeText = char(eventTypeNode.item(0).getTextContent);
                            if i~= 0
                                tmp = strsplit(eventTypeText, '|');
                                eventTypeText = tmp{2};
                                typeTooltip = tmp{1};
                            end
                            obj.EventConcepts{end+1} = eventConceptText;
                            if strcmp(eventTypeText, 'Stages') == 1
                                obj.EventStages{end+1} = eventConceptText;
                                SleepStageNames{end+1} = eventConceptText;
                            end
                            
                            %yh add spindle annotation 10/12/17 uncomment if spindle if an event type
                            %if strcmp(eventTypeText, 'Spindle') == 1
                            %    obj.EventStages{end+1} = eventConceptText;
                            %    sleepSpindleNames{end+1} = eventConceptText;
                            %end
                            %end yh 10/12/17
                            
                            if i ~= 0 %&& strcmp(eventTypeText, 'Stages') == 0
                                try
                                    InputChName = char(events.item(i).getElementsByTagName('SignalLocation').item(0).getTextContent);
                                    % is InputChName is empty, how to deal with it:
                                catch
                                    % No input tag found
                                end
                            end
                            
                            % Check if EventConcept contains these stages:
%                             Temp = strfind(eventConceptText,'Desaturation');
                            Temp = strfind(tooltip,'desaturation');
                            if ~isempty(Temp)
                                hasDesaturation = 1;
                                % change str2num to str2double todo
                                try
                                    nadirNode = events.item(i).getElementsByTagName('SpO2Nadir');
                                    nadirNum = str2num(nadirNode.item(0).getTextContent);
                                catch ex3
                                    obj = obj.logErr('Cannot found <SpO2Nadir> tag', i);
                                    continue
                                end
                                try
                                    baselineNode = events.item(i).getElementsByTagName('SpO2Baseline');
                                    baselineNum = str2num(baselineNode.item(0).getTextContent);
                                catch ex4
                                    obj = obj.logErr('Cannot found <SpO2Baseline> tag', i);
                                    continue
                                end
                            end
                        catch ex0
                            obj = obj.logErr('Cannot found <EventConcept> tag', i);                            
                            continue % ignore this event
                        end
                        try
                            startNode = events.item(i).getElementsByTagName('Start');
                            startNum = str2num(startNode.item(0).getTextContent);
                        catch ex1
                            obj = obj.logErr('Cannot found <Start> tag', i);
                            continue
                        end
                        try
                            durationNode = events.item(i).getElementsByTagName('Duration');
                            durationNum = str2num(durationNode.item(0).getTextContent);
                        catch ex2
                            obj = obj.logErr('Cannot found <Duration> tag', i);
                            continue
                        end
                        
                        if strcmp(eventTypeText, 'Stages') == 1
                            SleepStageDur{end+1} = durationNum;
                            SleepStageSt{end+1} = startNum;
                        end
                        
                        % Mod2
                        if strcmp(eventTypeText, 'Arousals') == 1
                            ArousalSt{end+1} = startNum;
                            ArousalDur{end+1} = durationNum;
                        end
                        
                        if strcmp(stagesNameVector(1),tooltip)==1
                            fprintf('sleep stage name: %s', eventConceptText);
                            SleepStageValues = [SleepStageValues, ones(1,durationNum)+3];
                        elseif strcmp(stagesNameVector(2),tooltip)==1
                            SleepStageValues = [SleepStageValues, ones(1,durationNum)+2];
                        elseif strcmp(stagesNameVector(3),tooltip)==1
                            SleepStageValues = [SleepStageValues, ones(1,durationNum)+1];
                        elseif strcmp(stagesNameVector(4),tooltip)==1
                            SleepStageValues = [SleepStageValues, ones(1,durationNum)];
                        elseif strcmp(stagesNameVector(5),tooltip)==1
                            SleepStageValues = [SleepStageValues, zeros(1,durationNum)];
                        elseif strcmp(stagesNameVector(6),tooltip)==1
                            SleepStageValues = [SleepStageValues, zeros(1,durationNum)+5];
                            % end
                        else
                        end
                        
                        %yh spindle annotation
                        if strcmp('Sleep Spindle',tooltip)==1
                            fprintf('sleep spindle name: %s', eventConceptText);
                            SleepSpindleValues = [SleepSpindleValues, zeros(1,startNum-length(SleepSpindleValues)), zeros(1,durationNum)+5.2];
                        end
                        %end yh
                        ithScoredEvent = struct(...
                            'EventType', eventTypeText, ...
                            'EventConcept', eventConceptText, ...
                            'Start', startNum, ...
                            'Duration', durationNum, ...
                            'InputCh', InputChName, ...
                            'SpO2Baseline', baselineNum, ...
                            'SpO2Nadir', nadirNum, ...
                            'ToolTip', tooltip ...
                        );
                        %------------------------------------- Validate Events
                        if i ~= 0
                            [eventValid, obj] = validateEvent(obj, ithScoredEvent, i);
                            if eventValid == 0
                                continue % error occurs during validation
                            end
                        end
                        %------------------------------ Construct Event Vector
                        if eventValid == 1
                            if i == 0
                                ScoredEvent(1) = ithScoredEvent;
                            else
                                % if this event is stage event, then do not
                                % include in the handles.ScoredEvent list
                                ScoredEvent(end+1) = ithScoredEvent;
                            end                            
                        end   
                        %------------------------------ End
                    end         
                    obj.EventConcepts = unique(obj.EventConcepts);
                    obj.EventStages = unique(obj.EventStages);
                end  % end if    
                fprintf('>>> Done.\n')
            end  % End Embedded Function                 
        end
        %--------------------------------------------- Validate Events
        function [isValid, obj] = validateEvent(obj, eventStruct, eventNum)
            %args: eventStruct is used for validating each of its fields
            %      eventNum is used for displaying error messages
            isValid = 1; % true
            eventErrMsg = '';
            if isempty(eventStruct.EventConcept)
                isValid = 0;
                eventErrMsg = strcat(eventErrMsg, 'Event name not found;');
            end
            if isempty(eventStruct.Start)
                isValid = 0;
                eventErrMsg = strcat(eventErrMsg, 'Start time empty;');
            end
            if isempty(eventStruct.Duration)
                isValid = 0;
                eventErrMsg = strcat(eventErrMsg, 'Duration empty;');
            end
            if strcmp(eventStruct.EventConcept, 'Desaturation')
                if isempty(eventStruct.SpO2Nadir)
                    isValid = 0;
                    eventErrMsg = strcat(eventErrMsg, 'SpO2Nadir empty;');
                end
                if isempty(eventStruct.SpO2Baseline)
                    isValid = 0;
                    eventErrMsg = strcat(eventErrMsg, 'SpO2Baseline empty;');
                end
            end
            if isValid == 0 && ~isempty(eventErrMsg)
                obj = obj.logErr(eventErrMsg, eventNum);
            end 
        end
        %-------------------------------------------- Log errors
        function obj = logErr(obj, message, eventNumber)
            %Error logging
              if isKey(obj.errMap, message)
                  obj.errMap(message) = [obj.errMap(message), ', ', num2str(eventNumber)];
              else 
                  obj.errMap(message) = [message, char(10), '  --> EventNumber: ', num2str(eventNumber)];
              end
        end
        
        function isStaging = isStageEvent(obj, scoredEvent)
            isStaging = 0;
            eventIndex = find(ismember(obj.EventStages, scoredEvent.EventConcept), 1);
            if eventIndex
                isStaging = 1;
            end
        end
    end
    %---------------------------------------------------- Private functions
    methods (Access=protected) 
    end
    %------------------------------------------------- Dependent Properties
    methods 
        %----------------------------------------------PhysioMiMi Variables
        %------------------------------------------------------- xmlEntries
        function value = get.xmlEntries(obj)
            value = obj.xmlEntriesP;
        end
        %------------------------------------------------------ ScoredEvent
        function value = get.ScoredEvent(obj)
            value = obj.ScoredEventP;
        end
        
        %-------------------------------------------------------- EventList
        function value = get.EventList(obj)
            value = obj.EventListP;
        end
        %------------------------------------------------------- EventTypes
        function value = get.EventTypes(obj)
            value = obj.EventTypesP;
        end
        %------------------------------------------------------- EventStart
        function value = get.EventStart(obj)
            value = obj.EventStartP;
        end                        
        %------------------------------------------------------ SleepStages
        function value = get.SleepStages(obj)
            value = obj.SleepStagesP;
        end
        function value = get.SleepStagesStart(obj)
            value = obj.SleepStagesSt;
        end
        function value = get.SleepStagesDuration(obj)
            value = obj.SleepStagesDur;
        end
        % mod2
        function value = get.ArousalsStart(obj)
            value = obj.ArousalSt;
        end
        function value = get.ArousalsDuration(obj)
            value = obj.ArousalDur;
        end
        %------------------------------------------------------ EpochLength
        function value = get.EpochLength(obj)
            value = obj.EpochLengthP;
        end
    end
    %------------------------------------------------- Dependent Properties
    methods(Static)        
        %---------------------------------------------------- GetEventTimes  
        function value = GetEventTimes(eventLabel, EventList, EventStart)
           % Return the time of the specified event
           
           % Define return value
           value = [];
           
           % Check for event typ
           eventIndex = strcmp(eventLabel, EventList);
           
           if ~isempty(eventIndex)
               value = EventStart(eventIndex);
           end 
        end
    end
end
% Staging,0,SRO:Wake,Epoch scored as Wake
% Staging,1,SRO:Stage1Sleep,Epoch scored as Stage 1 Sleep
% Staging,2,SRO:Stage2Sleep,Epoch scored as Stage 2 Sleep
% Staging,3,SRO:Stage3Sleep,Epoch scored as Stage 3 Sleep
% Staging,4,SRO:Stage4Sleep,Epoch scored as Stage 4 Sleep
% Staging,5,SRO:RapidEyeMovement,Epoch scored as REM Sleep
% Staging,6,SRO:MovementTime,Epoch scored as ppt time spent in movement
% Staging,9,SRO:UnscoredEpoch,Unscored epoch
% Staging,10,SRO:ArtifactEpoch,Epoch scored as artifact