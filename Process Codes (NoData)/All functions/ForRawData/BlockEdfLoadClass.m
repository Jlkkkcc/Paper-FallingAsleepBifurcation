classdef BlockEdfLoadClass
    %BlockEdfLoadClass Load multiple EDFs efficiently
    %   The class supports functional prototpes defined in the functional
    %   form (BlockEdfLoadClass). The class version of Block EDF loader is 
    %   designed to include additional functionality as compared to the 
    %   procedureal version. The new functionality is added in order to 
    %   reduce the amount of time required to analyze data.
    %   
    %   The EDF file name is the only required parameter. The user can 
    %   optionally set a signal list and a list of epochs. Setting a signal
    %   and epoch list can substranially reduce the amount of computer
    %   memory required to load/manipulated data stored in an EDF file.
    %
    %   Loaded data can either be stored in the class or passed directly to
    %   the calling function (without duplication), which allows for memory
    %   requirements to be managed by the user.
    %
    %   Functions for summarizing the header and signal header to the
    %   console are included.  A function for generating time series figures
    %   is included.
    %
    %   A procedure for identifying edf.* files recurrsively is included as
    %   a static function.
    % 
    %   A procedure for checking the EDF to the original specification is 
    %   provided.
    %
    %  The loader is designed to load the EDF file described in: 
    % 
    %    Bob Kemp, Alpo Värri, Agostinho C. Rosa, Kim D. Nielsen and John Gade 
    %    "A simple format for exchange of digitized polygraphic recordings" 
    %    Electroencephalography and Clinical Neurophysiology, 82 (1992): 
    %    391-393.%
    %
    %  An online description of the EDF format can be found at:
    %    http://www.edfplus.info/
    %
    %  Our EDF tools can be found at:
    %                  http://sleep.partners.org/edf/
    %
    %  Updated versions will be posted:
    %                  http://https://github.com/DennisDean/
    %
    %  Additional Documentation and examples will be posted:
    %                  http://sleepdata.org/tools
    %
    %  The header from the BlockEdfLoad file is copied below. Details
    %  specific to the class version follows.
    %
    % Public Properties:
    %     Required:
    %     edfFN : EDF file name with path
    %     
    %     Optional:
    %     signalLabels : Cell list of signal labels to load
    %     epochs:        [start epoch, end epoch] to load
    %     outArgClass:   1 to return class, 0 to return mimic return of 
    %                    BlockEdfLoad [header signal_header signal_cell]
    %     numCompToLoad: Sets the number of components to load. 
    %                         1. header
    %                         2. header and signal header
    %                         3. header, signal_header, and signal_cell
    %     tmax:          Display duration from start of signal
    %     fid:           Figure ids for figures created by class
    %     headerTxtFn:   Text file for saving edf header information
    %     signalHeaderTxtFn: 
    %                    Textfile for saving signal header information
    %     checkTxtFn:    Textfile for saving EDF check information
    %
    %  Dependent Properties:
    %     edf:           Structure holds loaded EDF components
    %     signal_labels: Returns the signal labels in the EDF file
    %     samples_in_record:
    %                    Number of samples per record
    %     sample_rate:   Sampling rate for each loaded signal
    %     tranducer_type: 
    %                     Transducer type contenst for each signal
    %     physical_dimension:
    %                     Physical dimensions contenst for each signal
    %     physical_max:   Physical maximums contenst for each signal
    %     physical_min:   Physical minimuns contenst for each signal
    %     digital_max:    Digital maximums contenst for each signal
    %     digital_min:    Digital minimuns contenst for each signal
    %     prefiltering:   Prefiltering contenst for each signal
    %     reserve_2:      Reserve 2 contenst for each signal
    %    
    %  Public Methods:
    %    Constructor:
    %       obj = BlockEdfLoadClass(edfFN)
    %       obj = BlockEdfLoadClass(edfFN, signalLabels) 
    %       obj = BlockEdfLoadClass(edfFN, signalLabels, epochs) 
    %    Load Prototypes (set load properties first)
    %       obj = obj.blockEdfLoad 
    %                 Default entire file, return class
    %       obj = obj.blockEdfLoad (outArgClass)
    %                 Select between class or structured return
    %       obj = obj.blockEdfLoad (outArgClass, numCompToLoad)
    %    Summary Functions
    %       obj.PrintEdfHeader
    %                 Write header contents to console
    %       obj.WriteEdfHeader
    %                 Write header to file defined in private properties
    %       obj.PrintEdfSignalHeader
    %                 Write signal header information to console
    %       obj.WriteEdfSignalHeader
    %                 Write signal header to file defined in private
    %                 properties
    %       obj.PlotEdfSignalStart
    %                 Create plot of initial signal for the intial duration
    %                 defined in the public properties (Default: 30
    %                 seconds)
    %       obj.CheckEdf
    %                 Check EDF header and signal header
    %       obj.DispCheck
    %                 Display results of check to console
    %       obj.WriteCheck
    %                 Write results of check to file
    %
    % Static Properties:
    %       GetEdfFileList:
    %                 Get list of EDF files in a folder (Multiple Forms)
    %       varargout = obj.GetEdfFileList
    %                 Default to current path
    %       varargout = obj.GetEdfFileList(folderPath)
    %                 Set folder path to search recursively
    %       obj.GetEdfFileList
    %                 Write file list to xls file defined internally
    %       fileList = obj.GetEdfFileList       
    %                 Cell array of EDF files
    %
    % Dependency:
    %       MS is required for option export of file list to an *.xls file
    %
    %
    % Acknowledgements
    %    Function uses DIRR function from MATLAB Central
    %    http://www.mathworks.com/matlabcentral/fileexchange/8682-dirr-find-files-recursively-filtering-name-date-or-bytes
    %
    %    Flatten file tree modeled after code available from MATLAB Central
    %
    %    EDF Checker modeled after EDFbrowser Code
    %    http://www.teuniz.net/edfbrowser/
    %
    %
    % ---------------------------------------------------------------------
    % blockEdfLoad Load EDF with memory block reads.
    % Function inputs an EDF file text string and returns the header,
    % header and each of the signals.
    %
    % Our EDF tools can be found at:
    %
    %                  http://sleep.partners.org/edf/
    %
    % The loader is designed to load the EDF file described in: 
    % 
    %    Bob Kemp, Alpo Värri, Agostinho C. Rosa, Kim D. Nielsen and John Gade 
    %    "A simple format for exchange of digitized polygraphic recordings" 
    %    Electroencephalography and Clinical Neurophysiology, 82 (1992): 
    %    391-393.
    %
    % An online description of the EDF format can be found at:
    % http://www.edfplus.info/
    %
    % Requirements:    Self contained, no external references 
    % MATLAB Version:  Requires R14 or newer, Tested with MATLAB 7.14.0.739
    %
    % Input (VARARGIN):
    %           edfFN : File text string 
    %    signalLabels : Cell array of signal labels to return (optional)
    %
    % Function Prototypes:
    %                                header = blockEdfLoad(edfFN)
    %                [header, signalHeader] = blockEdfLoad(edfFN)
    %    [header, signalHeader, signalCell] = blockEdfLoad(edfFN)
    %    [header, signalHeader, signalCell] = blockEdfLoad(edfFN, signalLabels)
    %    [header, signalHeader, signalCell] = blockEdfLoad(edfFN, signalLabels, epochs)
    %
    % Output (VARARGOUT):
    %          header : A structure containing variables for each header entry
    %    signalHeader : A structured array containing signal information, 
    %                   for each structure present in the data
    %      signalCell : A cell array that contains the data for each signal
    %
    % Output Structures:
    %    header:
    %       edf_ver
    %       patient_id
    %       local_rec_id
    %       recording_startdate
    %       recording_starttime
    %       num_header_bytes
    %       reserve_1
    %       num_data_records
    %       data_record_duration
    %       num_signals
    %    signalHeader (structured array with entry for each signal):
    %       signal_labels
    %       tranducer_type
    %       physical_dimension
    %       physical_min
    %       physical_max
    %       digital_min
    %       digital_max
    %       prefiltering
    %       samples_in_record
    %       reserve_2
    %
    % BlockEdfLoad Version: 0.1.16
    %
    %----------------------------------------------------------------------
    % BlockEdfLoadClass Description
    %
    %  Public Properties:
    %                      edfFN                 
    %               signalLabels  
    %                     epochs     
    %                outArgClass     
    %              numCompToLoad
    %            Display Options
    %                       tmax  
    %                        fid
    %
    %  Dependent Properties:
    %    
    %                        edf   EDF Header
    %              signal_labels   Signal Labels (cell array)
    %          samples_in_record   Samples in record 
    %                sample_rate   Sampling rate
    % 
    %          signalDurationSec     
    %  returnedSignalDurationSec     
    %             num30SecEpochs      
    %     returnedNum30SecEpochs  
    %              num5MinEpochs 
    %
    % Function Prototype
    %  Public:
    %     obj = BlockEdfLoadClass
    %     varargout = blockEdfLoad(obj, varargin)
    %     PrintEdfHeader(obj)
    %     PrintEdfHeader(obj)
    %     WriteEdfHeader(obj)
    %     WriteEdfHeader(obj)
    %     PrintEdfSignalHeader(obj)
    %     WriteEdfSignalHeader(obj)
    %     obj = PlotEdfSignalStart(obj)
    %
    % Static Functions
    %     varargout = GetEdfFileList(folder)
    %     function varargout = GetEdfFileListInfo(varargin)
    %
    % Version: 0.1.25
    %
    % ---------------------------------------------
    % Dennis A. Dean, II, Ph.D
    %
    % Program for Sleep and Cardiovascular Medicine
    % Brigam and Women's Hospital
    % Harvard Medical School
    % 221 Longwood Ave
    % Boston, MA  02149
    %
    % File created: April 29, 2013
    % Last updated: April 24, 2014 
    %    
    % Copyright © [2014] The Brigham and Women's Hospital, Inc. THE BRIGHAM AND 
    % WOMEN'S HOSPITAL, INC. AND ITS AGENTS RETAIN ALL RIGHTS TO THIS SOFTWARE 
    % AND ARE MAKING THE SOFTWARE AVAILABLE ONLY FOR SCIENTIFIC RESEARCH 
    % PURPOSES. THE SOFTWARE SHALL NOT BE USED FOR ANY OTHER PURPOSES, AND IS
    % BEING MADE AVAILABLE WITHOUT WARRANTY OF ANY KIND, EXPRESSED OR IMPLIED, 
    % INCLUDING BUT NOT LIMITED TO IMPLIED WARRANTIES OF MERCHANTABILITY AND 
    % FITNESS FOR A PARTICULAR PURPOSE. THE BRIGHAM AND WOMEN'S HOSPITAL, INC. 
    % AND ITS AGENTS SHALL NOT BE LIABLE FOR ANY CLAIMS, LIABILITIES, OR LOSSES 
    % RELATING TO OR ARISING FROM ANY USE OF THIS SOFTWARE.
    %
    
    %---------------------------------------------------- Public Properties
    properties (Access = public)
        % Required input
        edfFN                    % EDF file name
        
        % Optional Input
        signalLabels  = {};        % Cell list of signal labels to load
        epochs        = [];        % Start and end epoch to load
        
        % Load/return options
        outArgClass   = 1;         % Flag to determine class or data output
        numCompToLoad = 3;         % Load entire file
        
        % Display Options
        tmax = 30;                 % Display duration from start of signal
        fid = [];
        
        % Output file name
        headerTxtFn = 'header.txt'
        signalHeaderTxtFn = 'signalHeader.txt'   
        checkTxtFn = 'edfCheck.txt';
        
        % Error message summary
        errSummary = '';
        % Errors
        errList = [];
        
        % Operation parameters
        SWAP_MIN_MAX = 0;          % Swap if digMin > difMax; phyMin > phyMax
        INVERT_SWAP_MIN_MAX = 0;   % Invert signal if min-max are swapped
                                   % 1: Swapped, 0:Not Swapped
    end
     %------------------------------------------------ Dependent Properties
    properties (Dependent = true)
        % Access to EDF components
        edf                       % EDF Header
        signal_labels             % Signal Labels (cell array)
        samples_in_record         % Samples in record 
        sample_rate               % Sampling rate
        tranducer_type            % Transducer type
        physical_dimension        % Physical dimensions
        physical_max              % Physical maximums
        physical_min              % Physical minimuns
        digital_max               % Digital maximums
        digital_min               % Digital minimuns
        prefiltering              % Prefiltering
        reserve_2                 % Reserve 2
        
        % Error Checking Information
        errMsg                    % Error Mesasges from check
        mostSeriousErrValue       % Error Value      
        mostSeriousErrMsg         % Error Message
        totNumDeviations          % total Number of deviations
        deviationByType           % Number of deviations by class
        errSummaryMessages        % Error descrptions
        errSummaryLabel           % Error labels for tables
        
        % Properties supporting analysis
        signalDurationSec         % signal durations in seconds
        returnedSignalDurationSec % returned signal duration in seconds       
        num30SecEpochs            % number of complete 30 sec epochs in EDF
        returnedNum30SecEpochs    % num of returned 30 sec epochs (complete)
        num5MinEpochs             % number of five minute epochs
        returnedNum5MinEpochs     % # of returned five minute epochs (Complete)
    end
    %--------------------------------------------------- Private Properties
    properties (Access = protected)
        num_argsin                % number of arguments sent by user
        
        % Properties for saving loaded edf components
        header = [];                   % EDF Header
        signalHeader = [];             % Signal Header
        signalCell = [];               % Signal Cell Array
        
        % Check variables
        errMsgP = {};                  % Error messages from tests
        errMsgSerP = [];               % Error test codes
        
        % Error Types
        noErrorP = 1;
        inputError = 2;
        edfDeviationP = 3;
        scalingErrorP = 4;
        accessErrorP = 5;
        scalingAndDataAccessErrorP = 6;
        numericErrorTypesP = [1:6];
            
        % Error Messages
        errSummaryMessagesP = {...
            'No Error: File meets published EDF standard'; ...
            'BlockEdfLoadClass: Input Error';...
            'EDF Deviation: Applications may run'; ...
            % 'Scaling Error: May cause applications to fail'; ...
            'Scaling Warning: Applications may run'; ...
            'Fatal Error: Access Error'; ...
            'Fatal Error: Data Access Error and Scaling Error'  ...
        };
        errSummaryLabelP = {...
            'No Error'; ...
            'Input Error';...
            'Deviation Error'; ...
            'Scaling Error'; ...
            'Fatal Error'; ...
            'Fatal+Scaling Error'  ...
        };
    
        % Define check parameters
        MAX_SIGNALS = 512;
        DIG_MIN = -62768;       % digital min check: (-32768, -62768)
        DIG_MAX = 62767;        % digital max check  (32767,  62767)
        
    end
    %------------------------------------------------------- Public Methods
    methods
        %------------------------------------------------------ Constructor
        function obj = BlockEdfLoadClass(varargin)
            % Record number of arguments
            obj.num_argsin = nargin;
            
            % Process input
            if nargin == 1
                obj.edfFN = varargin{1};
            elseif nargin == 2
               obj.edfFN = varargin{1};
               obj.signalLabels = varargin{2};
            elseif nargin == 3 
               obj.edfFN = varargin{1};
               obj.signalLabels = varargin{2};  
               obj.epochs = varargin{3};
            else
                % Echo supported function prototypes to console
                fprintf('belObj = blockEdfLoad(edfFN)\n');
                fprintf('belObj = blockEdfLoad(edfFN, signalLabels)\n');
                fprintf('belObj = = blockEdfLoad(edfFN, signalLabels, epochs)\n');

                % Call MATLAB error function
                error('Function prototype not valid');
            end    
        end
        %--------------------------------------------------- Block EDF Load
        function varargout = blockEdfLoad(obj, varargin)
            % blockEdfLoad Load EDF with memory block reads.            

            %------------------------------------------------------------ Process input
            % Get load parameters
            edfFN = obj.edfFN;              % EDF files name
            signalLabels = obj.signalLabels;% Labels of signals to return
            epochs = obj.epochs;            % Start and end epoch to return

            % Get default load and return options
            outArgClass = obj.outArgClass;
            numCompToLoad = obj.numCompToLoad;
            
            % Process input
            if nargin == 2
                numCompToLoad = varargin{1};
            elseif nargin == 3;
                numCompToLoad = varargin{1};
                outArgClass = varargin{2};
            end

            % Set return value defaults
            header = [];
            signalHeader = [];
            signalCell = [];
            
            %-------------------------------------------------- Input check
            % Check that first argument is a string
            if   ~ischar(edfFN)
                msg = ('First argument is not a string.');
                obj.errMsgP(end+1) = {msg};
                error(msg);
            end
            % Check that first argument is a string
            if  ~iscellstr(signalLabels)
                msg = ('Second argument is not a valid text string.');
                obj.errMsgP(end+1) = {msg};
                error(msg);
            end
            % Check that first argument is a string
            if  and(obj.num_argsin ==3, length(epochs)~=2)
                msg = ('Specify epochs = [Start_Epoch End_Epoch.');
                obj.errMsgP(end+1) = {msg};
                error(msg);
            end

            %---------------------------------------  Load File Information
            % Load edf header to memory
            [fid, msg] = fopen(edfFN);

            % Proceed if file is valid
            if fid < 0
                % file id is not valid
                msg = sprintf('Could not open file: %s\n', edfFN);
                obj.errMsgP(end+1) = {msg};
                error(msg);    
            end


            % Open file for reading
            % Load file information not used in this version but will be used in
            % class version
            [filename, permission, machineformat, encoding] = fopen(fid);

            %-------------------------------------------------- Load Header
            try
                % Load header information in one call
                edfHeaderSize = 256;
                [A count] = fread(fid, edfHeaderSize);
            catch exception
                msg = 'Could not access header. Check available memory.';
                obj.errMsgP(end+1) = {msg};
                error(msg);
            end

            %----------------------------------------------------- Process Header Block
            % Create array/cells to create struct with loop
            headerVariables = {...
                'edf_ver';            'patient_id';         'local_rec_id'; ...
                'recording_startdate';'recording_starttime';'num_header_bytes'; ...
                'reserve_1';          'num_data_records';   'data_record_duration';...
                'num_signals'};
            headerVariablesConF = {...
                @strtrim;   @strtrim;   @strtrim; ...
                @strtrim;   @strtrim;   @str2num; ...
                @strtrim;   @str2num;   @str2num;...
                @str2num};
            headerVariableSize = [ 8; 80; 80; 8; 8; 8; 44; 8; 8; 4];
            headerVarLoc = vertcat([0],cumsum(headerVariableSize));
            headerSize = sum(headerVariableSize);

            % Create Header Structure
            header = struct();
            for h = 1:length(headerVariables)
                conF = headerVariablesConF{h};
                value = conF(char((A(headerVarLoc(h)+1:headerVarLoc(h+1)))'));
                header = setfield(header, headerVariables{h}, value);
            end
            
            num_data_records = header.num_data_records;
            

            % End Header Load section

            %------------------------------------------------------- Load Signal Header
            if numCompToLoad >= 2
                try 
                    % Load signal header into memory in one load
                    edfSignalHeaderSize = header.num_header_bytes - headerSize;
                    [A count] = fread(fid, edfSignalHeaderSize);
                catch exception
                    msg = 'Could not access signal header. Check available memory.';
                    obj.errMsgP(end+1) = {msg};
                    error(msg);
                end

                %------------------------------------------ Process Signal Header Block
                % Create arrau/cells to create struct with loop
                signalHeaderVar = {...
                    'signal_labels'; 'tranducer_type'; 'physical_dimension'; ...
                    'physical_min'; 'physical_max'; 'digital_min'; ...
                    'digital_max'; 'prefiltering'; 'samples_in_record'; ...
                    'reserve_2' };
                signalHeaderVarConvF = {...
                    @strtrim; @strtrim; @strtrim; ... 
                    @str2num; @str2num; @str2num; ...
                    @str2num; @strtrim; @str2num; ...
                    @strtrim };
                num_signal_header_vars = length(signalHeaderVar);
                num_signals = header.num_signals;
                signalHeaderVarSize = [16; 80; 8; 8; 8; 8; 8; 80; 8; 32];
                signalHeaderBlockSize = sum(signalHeaderVarSize)*num_signals;
                signalHeaderVarLoc = vertcat([0],cumsum(signalHeaderVarSize*num_signals));
                signalHeaderRecordSize = sum(signalHeaderVarSize);

                % Create Signal Header Struct
                signalHeader = struct(...
                    'signal_labels', {},'tranducer_type', {},'physical_dimension', {}, ...
                    'physical_min', {},'physical_max', {},'digital_min', {},...
                    'digital_max', {},'prefiltering', {},'samples_in_record', {},...
                    'reserve_2', {});

                % Get each signal header varaible
                for v = 1:num_signal_header_vars
                    varBlock = A(signalHeaderVarLoc(v)+1:signalHeaderVarLoc(v+1))';
                    varSize = signalHeaderVarSize(v);
                    conF = signalHeaderVarConvF{v};
                    for s = 1:num_signals
                        varStart = varSize*(s-1)+1;
                        varEnd = varSize*s;
                        value = conF(char(varBlock(varStart:varEnd)));

                        structCmd = ...
                            sprintf('signalHeader(%.0f).%s = value;',s, signalHeaderVar{v});
                        eval(structCmd);
                    end
                end
            end % End Signal Load Section

            %-------------------------------------------------------- Load Signal Block
            if numCompToLoad >= 3
                % Read digital values to the end of the file
                try
                    % Set default error mesage
                    errMsg = 'File load error. Check available memory.';

                    % Load strategy is dependent on input
                    if obj.num_argsin == 1
                        % Load entire file
                        [A count] = fread(fid, 'int16');
                    else 
                        % Get signal label information
                        edfSignalLabels = arrayfun(...
                            @(x)signalHeader(x).signal_labels, [1:header.num_signals],...
                                'UniformOutput', false);
                        signalIndexes = arrayfun(...
                            @(x)find(strcmp(x,edfSignalLabels)), signalLabels,...
                                'UniformOutput', false);

                        % Check that specified signals are present
                        signalIndexesCheck = cellfun(...
                            @(x)~isempty(x), signalIndexes, 'UniformOutput', false);
                        signalIndexesCheck = int16(cell2mat(signalIndexesCheck));
                        if sum(signalIndexesCheck) == length(signalIndexes);
                            % Indices are specified
                            signalIndexes = cell2mat(signalIndexes);
                        else
                            % Couldn't find at least one signal label
                            errMsg = 'Could not identify signal label';
                            error(errMsg);
                        end

                        edfSignalSizes = arrayfun(...
                            @(x)signalHeader(x).samples_in_record, [1:header.num_signals]);
                        edfRecordSize = sum(edfSignalSizes);

                        % Identify memory locations to record
                        endLocs = cumsum(edfSignalSizes)';
                        startLocs = [1;endLocs(1:end-1)+1];
                        signalLocs = [];
                        for s = signalIndexes
                            signalLocs = [signalLocs; [startLocs(s):1:endLocs(s)]'];
                        end
                        sizeSignalLocs = length(signalLocs);

                        % Load only required signals reduce memory calls
                        loadedSignalMemory = header.num_data_records*...
                            sum(edfSignalSizes(signalIndexes));
                        A = zeros(loadedSignalMemory,1);
                        for r = 1:header.num_data_records
                            [a count] = fread(fid, edfRecordSize, 'int16');
                            A([1+sizeSignalLocs*(r-1):sizeSignalLocs*r]) = a(signalLocs);
                        end

                        % Reset global varaibles, which enable reshape functions to
                        % work correctly
                        header.num_signals = length(signalLabels);
                        signalHeader = signalHeader(signalIndexes);
                        num_signals = length(signalIndexes);
                    end

                    %num_data_records
                catch exception
                    obj.errMsgP(end+1) = {errMsg};
                    error(errMsg);
                end
                %------------------------------------------------- Process Signal Block
                % Get values to reshape block
                num_data_records = header.num_data_records;
                getSignalSamplesF = @(x)signalHeader(x).samples_in_record;
                signalSamplesPerRecord = arrayfun(getSignalSamplesF,[1:num_signals]);
                recordWidth = sum(signalSamplesPerRecord);

                % Reshape - Each row is a data record
                A = reshape(A, recordWidth, num_data_records)';

                % Create raw signal cell array
                signalCell = cell(1,num_signals);
                signalLocPerRow = horzcat([0],cumsum(signalSamplesPerRecord));
                for s = 1:num_signals
                    % Get signal location
                    signalRowWidth = signalSamplesPerRecord(s);
                    signalRowStart = signalLocPerRow(s)+1;
                    signaRowEnd = signalLocPerRow(s+1);

                    % Create Signal
                    signal = reshape(A(:,signalRowStart:signaRowEnd)',...
                        signalRowWidth*num_data_records, 1);

                    % Get scaling factors
                    dig_min = signalHeader(s).digital_min;
                    dig_max = signalHeader(s).digital_max;
                    phy_min = signalHeader(s).physical_min;
                    phy_max = signalHeader(s).physical_max;

                    % Swap Digital Min-Max
                    if and((obj.SWAP_MIN_MAX == 1), dig_max < dig_min) 
                        signalHeader(s).digital_min = dig_max;
                        signalHeader(s).digital_max = dig_min;
                        dig_min = signalHeader(s).digital_min;
                        dig_max = signalHeader(s).digital_max;
                    end

                    % Swap Physical Min-Max
                    factor = 1; 
                    if and((obj.SWAP_MIN_MAX == 1), phy_max < phy_min) 
                        signalHeader(s).physical_min = phy_max;
                        signalHeader(s).physical_max = phy_min;
                        phy_min = signalHeader(s).physical_min;
                        phy_max = signalHeader(s).physical_max;
                        
                        % Invert signal if requested
                        if obj.INVERT_SWAP_MIN_MAX == 1
                            factor = -1;
                        end
                    end

                    
                    % Convert to analog signal
                    value = (signal-dig_min)/(dig_max-dig_min);
                    value = value.*double(phy_max-phy_min)+phy_min; 
                    
                    signalCell{s} = value*factor;
                end

            end
            
            %------------------------------------ Reduce signal if required
            % End Signal Load Section
            % Check if a reduce signal set is requested
            if ~isempty(epochs)
               % Determine signal sampling rate      
               signalSamples = arrayfun(...
                   @(x)signalHeader(x).samples_in_record, [1:num_signals]);
               signalIndex = ones(num_signals, 1)*[epochs(1)-1 epochs(2)]*30;
               samplesPerSecond = (signalSamples/header.data_record_duration)';
               signalIndex = signalIndex .* [samplesPerSecond samplesPerSecond];
               signalIndex(:,1) = signalIndex(:,1)+1;

               % Redefine signals to include specified epochs 
               signalIndex = int64(signalIndex);
               for s = 1:num_signals
                   signal = signalCell{s};
                   index = [signalIndex(s,1):signalIndex(s,2)];
                   signalCell{s} = signal(index);
               end
            end
               
            %------------------------------------------ Create return value
            % Check if object return requested
            if outArgClass == 1
               % Record edf information
               obj.header = header;
               obj.signalHeader = signalHeader;
               obj.signalCell = signalCell;
               
               % Assign output
               varargout{1} = obj;           
            else
                % Mirror functional form return
                if numCompToLoad < 2
                   varargout{1} = header;
                elseif numCompToLoad == 2
                   varargout{1} = header;
                   varargout{2} = signalHeader;
                elseif numCompToLoad == 3
                   % Create Output Structure
                   varargout{1} = header;
                   varargout{2} = signalHeader;
                   varargout{3} = signalCell;
                end % End Return Value Function       
            end
            
            % Close file explicitly
            if fid > 0 
                fclose(fid);
            end

        end % End of blockEdfLoad function
        %------------------------------------------------- Console Printing
        function PrintEdfHeader(obj)
            % Write header information to screen
            fprintf('Header:\n');
            fprintf('%30s:  %s\n', 'edf_ver', obj.header.edf_ver);
            fprintf('%30s:  %s\n', 'patient_id', obj.header.patient_id);
            fprintf('%30s:  %s\n', 'local_rec_id', obj.header.local_rec_id);
            fprintf('%30s:  %s\n', ...
                'recording_startdate', obj.header.recording_startdate);
            fprintf('%30s:  %s\n', ...
                'recording_starttime', obj.header.recording_starttime);
            fprintf('%30s:  %.0f\n', 'num_header_bytes', ...
                obj.header.num_header_bytes);
            fprintf('%30s:  %s\n', 'reserve_1', obj.header.reserve_1);
            fprintf('%30s:  %.0f\n', 'num_data_records', ...
                obj.header.num_data_records);
            fprintf('%30s:  %.0f\n', ...
                'data_record_duration', obj.header.data_record_duration);
            fprintf('%30s:  %.0f\n', 'num_signals', obj.header.num_signals);    
        end
        %--------------------------------------------------- Write Printing
        function WriteEdfHeader(obj)
            % Write header information to fil
            writeToFile = 1;
            if 1 == writeToFile
                fid = fopen(obj.headerTxtFn, 'w');
                
                fprintf(fid,'\n\nFile Name: %s\n%s):',obj.edfFN);  
                fprintf(fid, 'Header:\n');
                fprintf(fid, '%30s:  %s\n', 'edf_ver', obj.header.edf_ver);
                fprintf(fid, '%30s:  %s\n', 'patient_id', obj.header.patient_id);
                fprintf(fid, '%30s:  %s\n', 'local_rec_id', obj.header.local_rec_id);
                fprintf(fid, '%30s:  %s\n', ...
                    'recording_startdate', obj.header.recording_startdate);
                fprintf(fid, '%30s:  %s\n', ...
                    'recording_starttime', obj.header.recording_starttime);
                fprintf(fid, '%30s:  %.0f\n', 'num_header_bytes', ...
                    obj.header.num_header_bytes);
                fprintf(fid, '%30s:  %s\n', 'reserve_1', obj.header.reserve_1);
                fprintf(fid, '%30s:  %.0f\n', 'num_data_records', ...
                    obj.header.num_data_records);
                fprintf(fid, '%30s:  %.0f\n', ...
                    'data_record_duration', obj.header.data_record_duration);
                fprintf(fid, '%30s:  %.0f\n', 'num_signals', obj.header.num_signals);                   
                
            end 
                      
            % Close file
            fclose(fid);
        end
        %------------------------------------ Function PrintEdfSignalHeader
        function PrintEdfSignalHeader(obj)
            % Write signalHeader information to screen
            % Write signalHeader information to screen
            fprintf('\n\nSignal Header:');

            % Plot each signal
            for s = 1:obj.header.num_signals
                % Write summary for each signal
                fprintf('\n\n%30s:  %s\n', ...
                    'signal_labels', obj.signalHeader(s).signal_labels);
                fprintf('%30s:  %s\n', ...
                    'tranducer_type', obj.signalHeader(s).tranducer_type);
                fprintf('%30s:  %s\n', ...
                    'physical_dimension', ...
                    obj.signalHeader(s).physical_dimension);
                fprintf('%30s:  %.3f\n', ...
                    'physical_min', obj.signalHeader(s).physical_min);
                fprintf('%30s:  %.3f\n', ...
                    'physical_max', obj.signalHeader(s).physical_max);
                fprintf('%30s:  %.0f\n', ...
                    'digital_min', obj.signalHeader(s).digital_min);
                fprintf('%30s:  %.0f\n', ...
                    'digital_max', obj.signalHeader(s).digital_max);
                fprintf('%30s:  %s\n', ...
                    'prefiltering', obj.signalHeader(s).prefiltering);
                fprintf('%30s:  %.0f\n', ...
                    'samples_in_record', ...
                    obj.signalHeader(s).samples_in_record);
                fprintf('%30s:  %s\n', 'reserve_2', ...
                    obj.signalHeader(s).reserve_2);
            end
        end
        %------------------------------------ Function WriteEdfSignalHeader
        function WriteEdfSignalHeader(obj)           
            % Write file to header
            writeToFile = 1;
            if 1 == writeToFile
                fid = fopen(obj.signalHeaderTxtFn, 'w');
                
                fprintf(fid,'\n\nFile Name: %s\n%s):',obj.edfFN);   
                fprintf(fid,'Signal Header:\n');   
                
                % Plot each signal
                for s = 1:obj.header.num_signals
                    % Write summary for each signal
                    fprintf(fid,'\n\n%30s:  %s\n', ...
                        'signal_labels', obj.signalHeader(s).signal_labels);
                    fprintf(fid, '%30s:  %s\n', ...
                        'tranducer_type', obj.signalHeader(s).tranducer_type);
                    fprintf(fid, '%30s:  %s\n', ...
                        'physical_dimension', ...
                        obj.signalHeader(s).physical_dimension);
                    fprintf(fid, '%30s:  %.3f\n', ...
                        'physical_min', obj.signalHeader(s).physical_min);
                    fprintf(fid, '%30s:  %.3f\n', ...
                        'physical_max', obj.signalHeader(s).physical_max);
                    fprintf(fid, '%30s:  %.0f\n', ...
                        'digital_min', obj.signalHeader(s).digital_min);
                    fprintf(fid, '%30s:  %.0f\n', ...
                        'digital_max', obj.signalHeader(s).digital_max);
                    fprintf(fid, '%30s:  %s\n', ...
                        'prefiltering', obj.signalHeader(s).prefiltering);
                    fprintf(fid, '%30s:  %.0f\n', ...
                        'samples_in_record', ...
                        obj.signalHeader(s).samples_in_record);
                    fprintf(fid, '%30s:  %s\n', 'reserve_2', ...
                        obj.signalHeader(s).reserve_2);
                end
                
                % Close file
                fclose(fid);
            end
        end
        %-------------------------------------------- EDF Validate Function
        function obj = CheckEdf(obj)
            % Assign error types numeric and text designations
            noError = obj.noErrorP;
            edfDeviation = obj.edfDeviationP;
            scalingError = obj.scalingErrorP;
            accessError = obj.accessErrorP;
            scalingAndDataAccessError = obj.scalingAndDataAccessErrorP;
                    
            % Check EDF
            errorMSG = {};       
            errorSeverity = [];  

            % Check Header Variables
            if ~isempty(obj.header)  
                %------------------------------------- Process Header Block
                % Create array/cells to create struct with loop
                headerVariables = {...
                    'edf_ver';               'patient_id';...
                    'local_rec_id';          'recording_startdate';...
                    'recording_starttime';   'num_header_bytes'; ...
                    'reserve_1';             'num_data_records'; ...
                    'data_record_duration';  'num_signals'};
                headerVariablesConF = {...
                    @strtrim;   @strtrim;   @strtrim; ...
                    @strtrim;   @strtrim;   @str2num; ...
                    @strtrim;   @str2num;   @str2num;...
                    @str2num};
                headerVariableSize = [ 8; 80; 80; 8; 8; 8; 44; 8; 8; 4];
                
                % Variable types
                headerVariablesStr = {...
                    'edf_ver';                'patient_id';           ...           
                    'local_rec_id';           'recording_startdate';  ...
                    'recording_starttime';    'reserve_1' };              
                headerVariablesNum = {...
                    'num_header_bytes';       'num_data_records'; ...
                    'data_record_duration';   'num_signals'};   
                
                % Checking Function
                checkAsciiF = @(x) and(int8(x)<32, int8(x)>127);
                check_num_header_bytesF = @(x)x==256+256*obj.header.num_signals;
                check_num_data_recordsF = @(x)x==256;
                check_data_record_durationF = @(x)x==256;
                check_num_signalsF = @(x)and(x>0,x<=256); 
                check_periods = @(x)and(strcmp(x(3),'.'),strcmp(x(6),'.'));
                check_day = @(x)and(str2num(x(1:2))>=0,str2num(x(1:2))<=31);
                check_month = @(x)and(str2num(x(4:5))>=0,str2num(x(4:5))<=12);
                % check_year =
                % @(x)and(str2num(x(1:2))>0,str2num(x(7:8))<=31); %%%
                % 28.08.96 failed, changed on Jan, 29, 2015
                check_year = @(x)and(str2num(x(7:8))>=0,str2num(x(7:8))<=99);
                check_hour = @(x)and(str2num(x(1:2))>=0,str2num(x(1:2))<=23);
                check_minute = @(x)and(str2num(x(4:5))>=0,str2num(x(4:5))<=59);
                check_second = @(x)and(str2num(x(7:8))>=0,str2num(x(7:8))<=59);
                % check_second = @(x)and(str2num(x(1:2))>0,str2num(x(7:8))<=60);
                check_num_signal = @(x)and(x>0, x<=obj.MAX_SIGNALS);
                check_data_records = @(x)x>0;
                check_data_record_duration = @(x)x>0;
                check_digital_range = @(x)and(x>=obj.DIG_MIN, x<=obj.DIG_MAX);
                check_samples_in_record = @(x)x>=1;
                
                % Check ASCII outcome
                for t = 1:length(headerVariablesStr)
                   error = arrayfun(checkAsciiF, ...
                       getfield(obj.header, headerVariablesStr{t}));
                   if sum(error)> 0
                       errMsg = sprintf...
                           ('Header Error: Non ASCII character in ''%s''',...
                                headerVariablesStr{t}); 
                       errorMSG(end+1) = {errMsg};
                       errorSeverity(end+1) = edfDeviation;
                   end
                end
                
                % Check EDF Version
                if strcmp(obj.header.edf_ver,'0') == 0
                   errorMSG(end+1) = {'Header Error: Unknown version'};
                   errorSeverity(end+1) = edfDeviation;
                end
                
                % Check date strings
                dateChecks = ...
                    {check_periods, check_day, check_month, check_year};
                dateErrorMsg = { ...
                     'Header Error: Recording_startdate not standard';...
                     'Header Error: Recording startdate day is mis-specified';... 
                     'Header Error: Recording startdate month is mis-specified';...
                     'Header Error: Recording startdate year is mis-specified'};
                for c = 1:length(dateChecks)
                    check = dateChecks{c};
                    if check(obj.header.recording_startdate)==0
                        errorMSG(end+1) = {dateErrorMsg{c}};
                        errorSeverity(end+1) = edfDeviation;
                    end
                end
                            
                % Check start time
                dateChecks = ...
                    {check_periods, check_hour, check_minute, check_second};
                dateErrorMsg = { ...
                     'Header Error: Recording start time not standard';...
                     'Header Error: Recording start time hour is not standard';... 
                     'Header Error: Recording start time minute is not standard';...
                     'Header Error: Recording start time second is not standard'};
                for c = 1:length(dateChecks)
                    check = dateChecks{c};
                    if check(obj.header.recording_starttime)==0
                        errorMSG(end+1) = {dateErrorMsg{c}};
                        errorSeverity(end+1) = edfDeviation;
                    end
                end  
                
                % Check number of signals in header
                if check_num_signal(obj.header.num_signals)==0
                    errMsg = ...
                        sprintf('Header Error: Number of signals out of range (1-%.0f)',obj.MAX_SIGNALS);
                    errorMSG(end+1) = {errMsg};
                    errorSeverity(end+1) = accessError;
                end   
                
                % Check number of header bytes
                if check_num_header_bytesF(obj.header.num_header_bytes)==0
                    errMsg = ...
                        'Header Error: Number of header bytes should be 256+256*num_signals';
                    errorMSG(end+1) = {errMsg};
                    errorSeverity(end+1) = accessError;
                end   
                
                % Check data record
                if check_data_records(obj.header.num_data_records)==0
                    errMsg = ...
                        'Header Error: Number of data records > 0';
                    errorMSG(end+1) = {errMsg};
                    errorSeverity(end+1) = accessError;
                end  
                
                 % Check data record duration
                if check_data_record_duration(obj.header.data_record_duration)==0
                    errMsg = ...
                        'Header Error: Data record duration > 0';
                    errorMSG(end+1) = {errMsg};
                    errorSeverity(end+1) = accessError;
                end  
                
            end
                        
            % Store Error Message
            obj.errMsgP = errorMSG;
             
            % Check Signal Header Variables
            if and(~isempty(obj.signalHeader),obj.numCompToLoad>1)
                % Redefine ASCII function to conform to other checks
                checkAsciiF = @(x) and(int8(x)>=32, int8(x)<=126);
                % Check signal labels
                signal_labels = obj.signal_labels;
                for s = 1:obj.header.num_signals
                   error = arrayfun(checkAsciiF, signal_labels{s});
                   if sum(error)< length(signal_labels{s})
                       errMsg = sprintf...
                           ('Signal Header: Non ASCII character in signal, %s (signal %.0f)',...
                                signal_labels{s}, s); 
                       errorMSG(end+1) = {errMsg};
                       errorSeverity(end+1) = edfDeviation;
                   end
                end
                
                % Check transducer type
                tranducer_type = obj.tranducer_type;
                for s = 1:obj.header.num_signals
                   error = arrayfun(checkAsciiF, ...
                       getfield(obj.signalHeader(s), 'tranducer_type'));
                   if sum(error)< length(tranducer_type{s})
                       errMsg = sprintf...
                           ('Signal Header: Non ASCII character in transducer type, %s (signal %.0f)',...
                                signal_labels{s}, s); 
                       errorMSG(end+1) = {errMsg};
                       errorSeverity(end+1) = edfDeviation;
                   end
                end
              
                % Check prefiltering 
                prefiltering = obj.prefiltering;
                for s = 1:obj.header.num_signals
                   error = arrayfun(checkAsciiF, prefiltering{s});
                   if sum(error)< length(prefiltering{s})
                       errMsg = sprintf...
                           ('Signal Header: Non ASCII character in prefiltering field, %s (signal %.0f)',...
                                signal_labels{s}, s); 
                       errorMSG(end+1) = {errMsg};
                       errorSeverity(end+1) = edfDeviation;
                   end
                end
              
                % Check prefiltering 
                reserve_2 = obj.reserve_2;
                for s = 1:obj.header.num_signals
                   error = arrayfun(checkAsciiF, reserve_2{s});
                   if sum(error)< length(reserve_2{s})
                       errMsg = sprintf...
                           ('Signal Header: Non ASCII character in reserve_2 field, %s (signal %.0f)',...
                                signal_labels{s}, s); 
                       errorMSG(end+1) = {errMsg};
                       errorSeverity(end+1) = edfDeviation;
                   end
                end
                
                % Check physical dimensions
                physical_dimension = obj.physical_dimension;
                for s = 1:obj.header.num_signals
                   error = arrayfun(checkAsciiF, physical_dimension{s});
                   if sum(error)< length(physical_dimension{s})
                       errMsg = sprintf...
                           ('Signal Header: Non ASCII character in physical dimension, %s (signal %.0f)',...
                                signal_labels{s}, s); 
                       errorMSG(end+1) = {errMsg};
                       errorSeverity(end+1) = edfDeviation;
                   end
                end
 
                % Check physical min ~= physical max
                physical_max = obj.physical_max;
                physical_min = obj.physical_min;
                for s = 1:obj.header.num_signals
                   if physical_max(s)==physical_min(s)
                       errMsg = sprintf...
                           ('Signal Header: Physical min is equal to physical max, %s (signal %.0f)',...
                                signal_labels{s}, s); 
                       errorMSG(end+1) = {errMsg};
                       errorSeverity(end+1) = scalingError;
                   end
                end

                % Check physical min < physical max
                physical_max = obj.physical_max;
                physical_min = obj.physical_min;
                for s = 1:obj.header.num_signals
                   if ~(physical_max(s) > physical_min(s))
                       errMsg = sprintf...
                           ('Signal Header: Physical min is not less than physical max, %s (signal %.0f)',...
                                signal_labels{s}, s); 
                       errorMSG(end+1) = {errMsg};
                       errorSeverity(end+1) = scalingError;
                   end
                end
                
                % Check digital max range
                digital_max = obj.digital_max;
                for s = 1:obj.header.num_signals
                   if check_digital_range(digital_max(s))== 0
                       errMsg = sprintf...
                           ('Signal Header: Digital max is out of range, %s (signal %.0f)',...
                                signal_labels{s}, s); 
                       errorMSG(end+1) = {errMsg};
                       errorSeverity(end+1) = scalingError;
                   end
                end                
                
                % Check digital min range
                digital_min = obj.digital_min;
                for s = 1:obj.header.num_signals
                   if check_digital_range(digital_min(s))== 0
                       errMsg = sprintf...
                           ('Signal Header: Digital max is out of range, %s (signal %.0f)',...
                                signal_labels{s}, s); 
                       errorMSG(end+1) = {errMsg};
                       errorSeverity(end+1) = scalingError;
                   end
                end                       
                
                % Check digital min and max relationship
                digital_min = obj.digital_min;
                for s = 1:obj.header.num_signals
                   if digital_min(s) >= digital_max(s)
                       errMsg = sprintf...
                           ('Signal Header: Digital min is not less than digital max, %s (signal %.0f)',...
                                signal_labels{s}, s); 
                       errorMSG(end+1) = {errMsg};
                       errorSeverity(end+1) = scalingError;
                   end
                end
                
                % Check samples in records
                samples_in_record = obj.samples_in_record;
                for s = 1:obj.header.num_signals
                   if ~check_samples_in_record(samples_in_record(s))
                       errMsg = sprintf...
                           ('Signal Header: Samples in record is less than 1, %s (signal %.0f)',...
                                signal_labels{s}, s); 
                       errorMSG(end+1) = {errMsg};
                       errorSeverity(end+1) = accessError;
                   end
                end
                
            end
            
            % Store Error Message
            obj.errMsgP = errorMSG;
            obj.errList = errorMSG; %%% TODO
            obj.errMsgSerP = errorSeverity;
            
            if ~isempty(obj.errMsgSerP)
                obj.errSummary = obj.errSummaryMessages{max(obj.errMsgSerP)};
            end
            
            % Check Header Variables
            % No signal check performed
            if ~isempty(obj.signalCell)
                
            end            
        end
        %------------------------------------------------------ Write Check
        function DispCheck(obj)
            % Write error summary to console
            writeToFile = 1;
            if 1 == writeToFile
                fid = 1; % 1 redirect output to the console
                
                % Checking header 
                fprintf('Results of last check: %s\n', obj.edfFN);
                if isempty(obj.errMsg)
                    fprintf(fid,'\tNo error Messages found\n');
                else
                    for m = 1:length(obj.errMsg)
                        fprintf(fid,'\t%s\n',obj.errMsg{m});
                    end
                end    
                
                % Write Summary Statement
                fprintf(fid,'Most Serious Error From Last Check: %s\n', obj.edfFN);
                mostSeriousError = max(obj.errMsgSerP);
                if isempty(mostSeriousError)
                    fprintf(fid,'\t%s\n', 'No Errors Found!');
                else
                    errSummaryMessages = obj.errSummaryMessagesP;
                    fprintf(fid,'\t%s\n', errSummaryMessages{mostSeriousError});
                end
            end 
        end
        %------------------------------------------------------ Write Check
        function WriteCheck(obj)
            % Write header information to fil
            writeToFile = 1;
            if 1 == writeToFile
                fid = fopen(obj.checkTxtFn, 'w');
                
                % Checking header 
                fprintf(fid,'Results of last check: %s\n', obj.edfFN);
                if isempty(obj.errMsg)
                    fprintf(fid,'\tNo error Messages found\n');
                else
                    for m = 1:length(obj.errMsg)
                        fprintf(fid,'\t%s\n',obj.errMsg{m});
                    end
                end    
                
                % Write Summary Statement
                fprintf(fid,'\nSummary of last check: %s\n', obj.edfFN);
                mostSeriousError = max(obj.errMsgSerP);
                errSummaryMessages = obj.errSummaryMessagesP;
                fprintf(fid,'\t%s\n', errSummaryMessages{mostSeriousError});
            end 
                      
            % Close file
            fclose(fid);
        end
        %----------------------------------------------- Plotting Functions
        %----------------------------------------------- PlotEdfSignalStart   
        function obj = PlotEdfSignalStart(obj)
            % Function for plotting start of edf signals
            
            % Create figure
            fid = figure();

            % Get number of signals
            num_signals = obj.header.num_signals;

            % Add each signal to figure
            for s = 1:num_signals
                % get signal
                signal =  obj.signalCell{s};
                record_duration = obj.header.data_record_duration;
                samplingRate = obj.signalHeader(s).samples_in_record/...
                    record_duration;    
                t = [0:length(signal)-1]/samplingRate';
                

                % Identify first 30 seconds
                indexes = find(t<=obj.tmax);
                signal = signal(indexes);
                t = t(indexes);

                % Normalize signal
                sigMin = min(signal);
                sigMax = max(signal);
                signalRange = sigMax - sigMin;
                signal = (signal - sigMin);
                if signalRange~= 0
                    signal = signal/(sigMax-sigMin); 
                end
                signal = signal -0.5*mean(signal) + num_signals - s + 1;         

                % Plot signal
                plot(t(indexes), signal(indexes));
                hold on
            end

            % Set title
            title(obj.header.patient_id);

            % Set axis limits
            v = axis();
            v(1:2) = [0 obj.tmax];
            v(3:4) = [-0.5 num_signals+1.5];
            axis(v);

            % Set x axis
            xlabel('Time(sec)');

            % Set yaxis labels
            signalLabels = cell(1,num_signals);
            for s = 1:num_signals
                signalLabels{num_signals-s+1} = ...
                    obj.signalHeader(s).signal_labels;
            end
            set(gca, 'YTick', [1:1:num_signals]);
            set(gca,'YTickLabel', signalLabels);

            % Save figure id
            obj.fid = fid;
            
        end        
    end
    %---------------------------------------------------- Private functions
    methods (Access=protected)
    end
    %------------------------------------------------- Dependent Properties
    methods
        %-------------------------------------------------------------- edf
        function value = get.edf(obj)
            % returns loaded edf components in a single structure
            value.header = obj.header;
            if obj.numCompToLoad >=2
                value.signalHeader = obj.signalHeader;
            end
            if obj.numCompToLoad >= 3
                value.signalCell = obj.signalCell;
            end
        end
        %---------------------------------------------------- signal_labels
        function value = get.signal_labels(obj)
            % returns loaded edf components in a single structure
            value = {};
            if obj.numCompToLoad >=2
                N = obj.header.num_signals;
                signalHeader = obj.signalHeader;
                value = arrayfun(@(x)signalHeader(x).signal_labels,[1:N],...
                    'UniformOutput', false);
            end
        end
        %------------------------------------------------------ sample_rate
        function value = get.sample_rate(obj)
            % returns loaded edf components in a single structure
            value = {};
            if obj.numCompToLoad >=2
                N = obj.header.num_signals;
                signalHeader = obj.signalHeader;
                value = arrayfun(@(x)signalHeader(x).samples_in_record,...
                    [1:N])...
                        /obj.header.data_record_duration;
            end
        end
        %--------------------------------------------------- tranducer_type
        function value = get.tranducer_type(obj)
            % returns loaded edf components in a single structure
            value = {};
            if obj.numCompToLoad >=2
                N = obj.header.num_signals;
                signalHeader = obj.signalHeader;
                value = ...
                    arrayfun(@(x)signalHeader(x).tranducer_type,[1:N], ...
                         'UniformOutput', 0);
            end
        end
       %------------------------------------------------ physical_dimension
        function value = get.physical_dimension(obj)
            % returns loaded edf components in a single structure
            value = {};
            if obj.numCompToLoad >=2
                N = obj.header.num_signals;
                signalHeader = obj.signalHeader;
                value = ...
                    arrayfun(@(x)signalHeader(x).physical_dimension,[1:N], ...
                         'UniformOutput', 0);
            end
        end
       %------------------------------------------------------ physical_max
       function value = get.physical_max(obj)
            % returns loaded edf components in a single structure
            value = [];
            if obj.numCompToLoad >=2
                N = obj.header.num_signals;
                signalHeader = obj.signalHeader;
                value = ...
                    arrayfun(@(x)signalHeader(x).physical_max,[1:N]);
            end
       end
       %------------------------------------------------------ physical_min
       function value = get.physical_min(obj)
            % returns loaded edf components in a single structure
            value = [];
            if obj.numCompToLoad >=2
                N = obj.header.num_signals;
                signalHeader = obj.signalHeader;
                value = ...
                    arrayfun(@(x)signalHeader(x).physical_min,[1:N]);
            end
       end
       %------------------------------------------------------- digital_max
       function value = get.digital_max(obj)
            % returns loaded edf components in a single structure
            value = [];
            if obj.numCompToLoad >=2
                N = obj.header.num_signals;
                signalHeader = obj.signalHeader;
                value = ...
                    arrayfun(@(x)signalHeader(x).digital_max,[1:N]);
            end
       end
       %------------------------------------------------------- digital_min
       function value = get.digital_min(obj)
            % returns loaded edf components in a single structure
            value = [];
            if obj.numCompToLoad >=2
                N = obj.header.num_signals;
                signalHeader = obj.signalHeader;
                value = ...
                    arrayfun(@(x)signalHeader(x).digital_min,[1:N]);
            end
       end
       %------------------------------------------------------ prefiltering
       function value = get.prefiltering(obj)
            % returns loaded edf components in a single structure
            value = [];
            if obj.numCompToLoad >=2
                N = obj.header.num_signals;
                signalHeader = obj.signalHeader;
                value = ...
                    arrayfun(@(x)signalHeader(x).prefiltering,[1:N],...
                        'UniformOutput', 0);
            end
       end
       %--------------------------------------------------------- reserve_2
       function value = get.reserve_2(obj)
            % returns loaded edf components in a single structure
            value = [];
            if obj.numCompToLoad >=2
                N = obj.header.num_signals;
                signalHeader = obj.signalHeader;
                value = ...
                    arrayfun(@(x)signalHeader(x).reserve_2,[1:N],...
                        'UniformOutput', 0);
            end
       end
       %------------------------------------------------- samples_in_record
       function value = get.samples_in_record(obj)
            % returns loaded edf components in a single structure
            value = [];
            if obj.numCompToLoad >=2
                N = obj.header.num_signals;
                signalHeader = obj.signalHeader;
                value = ...
                    arrayfun(@(x)signalHeader(x).samples_in_record,[1:N],...
                        'UniformOutput', 1);
            end
       end   
       %----------------------------------------------------- Check Results
       %------------------------------------------------------------ errMsg
       function value = get.errMsg(obj)
            % returns loaded edf components in a single structure
            value = obj.errMsgP;
       end
       %------------------------------------------------------------ errMsg
       function value = get.mostSeriousErrValue(obj)
            % returns loaded edf components in a single structure
            value = max(obj.errMsgSerP);
       end
       %------------------------------------------------------------ errMsg
       function value = get.mostSeriousErrMsg(obj)
            % returns loaded edf components in a single structure
            if ~isempty(obj.errMsgSerP);
                value = obj.errSummaryMessagesP{max(obj.errMsgSerP)};
            else
                value = 'No error found!';
            end
       end
       %-------------------------------------------------- totNumDeviations
       function value = get.totNumDeviations(obj)
            % returns loaded edf components in a single structure
            if ~isempty(obj.errMsgSerP)
            	value = 0;
            else
                value = length(obj.errMsgSerP);
            end       
       end
       %--------------------------------------------------- deviationByType
       function value = get.deviationByType(obj)
            % returns loaded edf components in a single structure
            numericErrorTypesP = obj.numericErrorTypesP;
            errMsgSerP = obj.errMsgSerP;
            countErrTypeF = @(x)sum(errMsgSerP == x);
            value = arrayfun(countErrTypeF, numericErrorTypesP);
            if sum(value) == 0
                % Set no error index to 1
                value(1) = 1;
            end
       end       
       %------------------------------------------------ errSummaryMessages
       function value = get.errSummaryMessages(obj)
            % returns loaded edf components in a single structure
            value = obj.errSummaryMessagesP;
       end
       %--------------------------------------------------- errSummaryLabel
       function value = get.errSummaryLabel(obj)
            % returns loaded edf components in a single structure
            value = obj.errSummaryLabelP;
       end
       %----------------------------------- Properties Supporting Analysis
        %------------------------------------------------ signalDurationSec
        function value = get.signalDurationSec(obj)
            % returns loaded edf components in a single structure
            value = [];
            if obj.numCompToLoad >= 2
                num_data_records = obj.header.num_data_records;
                data_record_duration = obj.header.data_record_duration;
                
                value = data_record_duration*num_data_records;
            end
        end
        %---------------------------------------- returnedSignalDurationSec
        function value = get.returnedSignalDurationSec(obj)
            % returns loaded edf components in a single structure
            value = [];
            if isempty(obj.epochs)
                value = obj.signalDurationSec;
            else
                value = (obj.epochs(2)-obj.epochs(1)+1)*30;
            end
        end
         %-------------------------------------------------- num30SecEpochs
        function value = get.num30SecEpochs(obj)
             % number of complete 30 sec epochs in EDF
            value = [];
            if obj.numCompToLoad >= 2
                num_data_records = obj.header.num_data_records;
                data_record_duration = obj.header.data_record_duration;
                
                value = floor((data_record_duration*num_data_records)/30);
            end
        end
        %------------------------------------------- returnedNum30SecEpochs
        function value = get.returnedNum30SecEpochs(obj)
            % num of returned 30 sec epochs (complete)
            value = [];
            if isempty(obj.epochs)
                num_data_records = obj.header.num_data_records;
                data_record_duration = obj.header.data_record_duration;
                value = floor((data_record_duration*num_data_records)/30);
            else
                value = (obj.epochs(2)-obj.epochs(1)+1);
            end
        end
        %---------------------------------------------------- num5MinEpochs
        function value = get.num5MinEpochs(obj)
            % number of five minute epochs
            value = [];
            if obj.numCompToLoad >= 2
                num_data_records = obj.header.num_data_records;
                data_record_duration = obj.header.data_record_duration;
                value = floor((data_record_duration*num_data_records)/300);
            end
        end
        %-------------------------------------------- returnedNum5MinEpochs
        function value = get.returnedNum5MinEpochs(obj)
            % number of returned five minute epochs (Complete)
            value = [];
            if isempty(obj.epochs)
                num_data_records = obj.header.num_data_records;
                data_record_duration = obj.header.data_record_duration;
                value = floor((data_record_duration*num_data_records)/300);
            else
                value = floor((obj.epochs(2)-obj.epochs(1)+1)/10);
            end
        end
    end
    %------------------------------------------------- Dependent Properties
    methods(Static)
        %--------------------------------------------------- GetEdfFileList
        function varargout = GetEdfFileList(folder)
            % Get EDF file information
            fileListCellwLabels = ...
                BlockEdfLoadClass.GetEdfFileListInfo(folder);
            fileListCell = fileListCellwLabels(2:end,:);
            
            % Generate file list 
            fileList = arrayfun(...
                @(x)strcat(fileListCell{x,end},'\',fileListCell{x,1}), ...
                    [1:size(fileListCell,1)], 'UniformOutput', false);
            fn = fileListCell(:,1);
             
            % Return content determined by number of calling arguments    
            if nargout == 1
                varargout{1} = fileList;
            elseif nargout == 2
                varargout{1} = fileList;
                varargout{2} = fn;
            else
                fprintf('filelist = obj.GetEdfFileList(folder)\n');
                fprintf('[filelist fn] = obj.GetEdfFileList(folder)\n');
                msg = 'Number of output arguments not supported';
                error(msg);
            end
        end
        %----------------------------------------------- GetEdfFileListInfo
        function varargout = GetEdfFileListInfo(varargin)
            % Create default value
            value = [];
            folderPath = '';
            xlsOut = 'edfFileList.xls';
            
            % Process input
            if nargin ==0
                % Open window
                folderPath = uigetdir(cd,'Set EDF search folder');    
                if folderPath == 0
                    error('User did not select folder');
                end
            elseif nargin == 1
                % Set EDF search path
                folderPath = varargin{1};
            else
                fprintf('fileStruct = obj.locateEDFs(path| )\n');
            end

            % Get File List
            fileTree  = dirr(folderPath, '\.edf');
            [fileList fileLabels]= flattenFileTree(fileTree, folderPath);
            fileList = [fileLabels;fileList];
            
            % Write output to xls file
            if nargout == 0
                xlsOut = strcat(folderPath, '\', xlsOut);
                xlswrite('edfFileList.xls',[fileLabels;fileList]);
            else
                varargout{1} = fileList;
            end
            
            %---------------------------------------------- FlattenFileTree
            function varargout = flattenFileTree(fileTree, folder)
                % Process recursive structure created by dirr (See MATLAB Central)
                % find directory and file entries
                dirMask = arrayfun(@(x)isstruct(fileTree(x).isdir) == 1, ...
                    [1:length(fileTree)]);
                fileMask = ~dirMask;

                % Recurse on each directory entry
                fileListD = {};
                if sum(int16(dirMask)) > 0
                   dirIndex = find(dirMask);
                   for d = dirIndex
                       folderR = strcat(folder,'\',fileTree(d).name);
                       fileListR = flattenFileTree(fileTree(d).isdir, folderR);
                       fileListD = [fileListD; fileListR];
                   end 
                end

                % Merge current and recursive list
                fileList = {};
                if sum(int16(fileMask)) > 0
                   fileIndex = find(fileMask);
                   for f = fileIndex
                       entry = {fileTree(f).name ...
                                fileTree(f).date  ...
                                fileTree(f).bytes  ...
                                fileTree(f).datenum ...
                                folder};
                       fileList = [fileList; entry];
                   end   
                end

                % Merg diretory and file list
                fileList = [fileList; fileListD];

                % Pass file list labels on export
                if nargout == 1
                    varargout{1} = fileList;
                elseif nargout == 2
                    varargout{1} = fileList;
                    varargout{2} = ...
                        {'name', 'date', 'bytes',  'datenum', 'folder'};
                end
            end
        end
    end    
end