function [header signalHeader signalCell] = blockEdfLoad(edfFN)
% blockEdfLoad Load EDF with memory block reads.
% File reads an EDF file text string and returns the header, signal header,
% and each of the signals.
%
% The loader is designed to load the EDF file described in: 
% 
%    Bob Kemp, Alpo Värri, Agostinho C. Rosa, Kim D. Nielsen and John Gade 
%    "A simple format for exchange of digitized polygraphic recordings" 
%    Electroencephalography and Clinical Neurophysiology, 82 (1992): 
%    391-393.
%
% An online view of the EDF format can be found: http://www.edfplus.info/
%
% Requirements:    Self contained, no external references 
% MATLAB Version:  Uses abstract functions, tested in MATLAB 7.14.0.739
%
% Input:
%           edfFN : File text string 
%
% Output:
%          header : A structure containing varaibles for each header entry
%    signalHeader : A structured arrau containing signal information, 
%                   for each structure present in the data
%      signalCell : A cell array that contains the data for each signal
%
% Structures:
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
% Version: 0.1.0
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
% File created: October 23, 2012
% Last updated:  December 6, 2012 
%    
% Copyright © [2012] The Brigham and Women's Hospital, Inc. THE BRIGHAM AND 
% WOMEN'S HOSPITAL, INC. AND ITS AGENTS RETAIN ALL RIGHTS TO THIS SOFTWARE 
% AND ARE MAKING THE SOFTWARE AVAILABLE ONLY FOR SCIENTIFIC RESEARCH 
% PURPOSES. THE SOFTWARE SHALL NOT BE USED FOR ANY OTHER PURPOSES, AND IS
% BEING MADE AVAILABLE WITHOUT WARRANTY OF ANY KIND, EXPRESSED OR IMPLIED, 
% INCLUDING BUT NOT LIMITED TO IMPLIED WARRANTIES OF MERCHANTABILITY AND 
% FITNESS FOR A PARTICULAR PURPOSE. THE BRIGHAM AND WOMEN'S HOSPITAL, INC. 
% AND ITS AGENTS SHALL NOT BE LIABLE FOR ANY CLAIMS, LIABILITIES, OR LOSSES 
% RELATING TO OR ARISING FROM ANY USE OF THIS SOFTWARE.
%

%---------------------------------------------------  Load File Information
% Load edf header to memory
[fid, msg] = fopen(edfFN);

% Proceed if file is valid
if fid <0
    % file id is not valid
    error(msg);    
end

% Load file information
% Not used in this version but will be used in class version
[filename, permission, machineformat, encoding] = fopen(fid);

%-------------------------------------------------------------- Load Header
edfHeaderSize = 256;
[A count] = fread(fid, edfHeaderSize);

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

%------------------------------------------------------- Load Signal Header
% Load signal header into memory
edfSignalHeaderSize = header.num_header_bytes - headerSize;
[A count] = fread(fid, edfSignalHeaderSize);

%---------------------------------------------- Process Signal Header Block
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

%-------------------------------------------------------- Load Signal Block
% Read digital values to the end of the record
[A count] = fread(fid, 'int16');

%----------------------------------------------------- Process Signal Block
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
    
    % Convert to analog signal
    value = double(signal) - (dig_max+dig_min)/2;
    value = value./(dig_max-dig_min);
    if phy_min >0
        value = -value;
    end
    signalCell{s} = value;
end

end