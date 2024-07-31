function [xRRI,fsRRI]=ECG_to_RRI(varargin)
% -------------------------------------------------------------------------
% ECG_to_RRI.m
% Author: David Looney
% Email: daveloon@gmail.com
% Date: Jan 2015
% 
% -------------------------------------------------------------------------
% Syntax: 
%          [xRRI,fsRRI]=ECG_to_RRI(xECG,fsECG);
% 
% -------------------------------------------------------------------------
% Description:
% 
% Converts ECG time series input to RR interval time series.
% 
%   Step 1: The ECG is bandpass filtered (suggested range 5 - 20 Hz)
%   Step 2: The R peaks are estimated. R peaks are expected to be larger
%           than parameter 'ampthresh', and it is expected that a minimum 
%           time interval exists between successive peaks 'timethresh'.
%   Step 3: The user is presented with a plot illustrating peaks that have
%           been identified (upper subplot) and the time difference between 
%           successive peaks (lower subplot). 
%           The user is then prompted with: 
% 
%           Change parameters (Y/N) ?
% 
%           Entering 'Y' allows the user to change the parameters before 
%           proceeding to the next step. If only a small number of peaks 
%           have been misclassified, the user may review these errors in the
%           next step by entering 'N'.
% 
%   Step 4: If anomalies have been detected, the user is then prompted
%           with:
% 
%           Possible anomaly detected on or after time x, remove (Y/N)? 
% 
%           For each detected anomaly, the user can then manually review
%           the plots to determine whether the anomaly is genuine (a peak
%           has not been detected OR a spurious peak has been detected). If
%           the anomaly is genuine, the user should enter 'Y' and the
%           anomaly will be corrected. Otherwise, the user should enter 'N'
%           and no changes will be made concerning the potential anomaly.
%           The sensitivity of the anomaly detection algorithm is controlled
%           by 'anomalyparam'.
%   Step 5: The difference between R peaks is sampled at regular intervals
%           (controlled by parameter fsRRI) to generate the RRI.
%   Step 6: The user is presented with a plot of the final RRI with
%           detected anomalies, and a plot of the RRI with anomalies removed. 
% 
% -------------------------------------------------------------------------
% Required inputs:
% 
%       xECG  ~  An ECG time series. The expected voltage range is microvolt.
%                At least 5 s of ECG is required to estimate the RRI.
%       fsECG ~  The sampling frequency of the ECG. Recordings with
%                sampling frequencies below 32 Hz will be rejected.
% 
% -------------------------------------------------------------------------
% Outputs:
% 
%       xRRI  ~  The RR interval (RRI) time series, sampled at regular
%                intervals.
%       fsRRI ~  The sampling frequency of the RRI.
% 
% -------------------------------------------------------------------------
% Parameters:
% 
%       Hp             ~    High pass ECG filter cutoff (in Hz).
%       Lp             ~    Low pass ECG filter cutoff (in Hz).
%       ampthresh      ~    Minimum R peak amplitude.
%       ampthresh_autodetect ~ 'Y' (default) or 'N'. The user can let the function 
%                           pick an appropriate value for 'ampthresh'.   
%       timethresh     ~    Minimum time distance between R peaks (in s). 
%       fsRRI          ~    The sampling frequency of the RRI.
%       anomalyparam   ~    Controls the sensitivity of the anomaly detection
%                           algorithm. The larger the value, the less
%                           sensitive the algorithm.
%           
% 
% Parameter values can be changed from their default values in the following
% way:
%          [xRRI,fsRRI]=ECG_to_RRI(xECG,fsECG,'timethresh',0.5);
% 
% In the above, the user has selected the parameter 'timethresh' as 0.5 s.
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
FontSize = 16; 
% load inputs and parameters
[xECG,fsECG,hpECG,lpECG,fsRRI,ampthresh,timethresh,anomalyparam,ampthresh_autodetect] = init(varargin{:});

% perform inversion test
inversion_test(xECG,fsECG);

% generate filter coefficients
try
    [b,a]=butter(4,[hpECG,lpECG]/(fsECG/2),'bandpass');
catch
    % For more recent releases of MATLAB
    [b,a]=butter(4,[hpECG,lpECG]/(fsECG/2),'pass');
end
% filter ECG
yECG=filtfilt(b,a,xECG);

% if ampthresh_autodetect is selected, the amplitude threshold is selected
% in an automatic fashion
if strcmp(ampthresh_autodetect,'Y')
    ampthresh=3.4*std(yECG);
end

% calculate RRI
[~,xRRI]=get_RRI(yECG,fsECG,fsRRI,ampthresh,timethresh,anomalyparam);



function inversion_test(x,fs)

% generate highpass filter coefficients
[b,a]=butter(4,[1]/(fs/2),'high');

y=filtfilt(b,a,x);

if skewness(y)<0
    disp('Warning. The ECG may be inverted.')
end



function [TQ,RR]=get_RRI(ECG,ECG_fs,RRI_fs,ampthresh,timethresh,anomalyparam)

FontSize = 16; 

T=[1:length(ECG)]/ECG_fs;

TQ=[1:T(end)*RRI_fs]/RRI_fs;

SATISFIED_WITH_PARAMETERS=0;
ANOMALIES_DETECTED=0;

while ~SATISFIED_WITH_PARAMETERS
    % find peaks in ECG data
    [~,T_peak]=ECG_peak_detection_v2(ECG,T,ECG_fs,ampthresh,timethresh);
    ip_str = input(['Change parameters (Y/N) ? '],'s');
    ip_str=upper(ip_str);
    if strcmp(ip_str,'Y')
        param_str = input(['Which parameter (ampthresh,timethresh,sign) ? '],'s');
        
        if strcmp(param_str,'ampthresh')
            
            ampthresh=input(['Enter ampthresh ']);
            
        elseif strcmp(param_str,'timethresh')
            
            timethresh=input(['Enter timethresh ']);
            
        elseif strcmp(param_str,'sign')
            
            ECG=ECG*-1;
            disp('Sign of ECG has been changed.')
                  
        else
            disp('Input not recognized.');
        end
            
        
    elseif strcmp(ip_str,'N')
        SATISFIED_WITH_PARAMETERS=1;
    else 
        disp('Input not recognized.');
    end
end

RR_unfiltered=((diff((T_peak))));

NN=30;


Count=2;
RR_filtered(1)=RR_unfiltered(1);
T_peak_filtered(1)=T_peak(1);

for n=2:1:length(RR_unfiltered)

    RR_median=median(RR_filtered(Count-1:-1:max(1,Count-NN)));

    if (RR_unfiltered(n)>((1+anomalyparam)*RR_median))||(RR_unfiltered(n)<(max((1-anomalyparam),0)*RR_median))
        USER_SATISFACTION=0;
        ANOMALIES_DETECTED=1;
        while ~USER_SATISFACTION
            user_str = input(['Possible anomaly detected on or after time ',num2str(T_peak(n)),' , remove (Y/N)? '],'s');
            user_str=upper(user_str);
            if strcmp(user_str,'Y')
                disp('Anomaly corrected.')
                USER_SATISFACTION=1;
            elseif strcmp(user_str,'N')
                disp('No change made.')
                RR_filtered(end+1)=RR_unfiltered(n);
                T_peak_filtered(end+1)=T_peak(n);
                Count=Count+1;

                
                USER_SATISFACTION=1;
            else
                disp('Input not recognized.');
            end
        end             
                    
    else

         RR_filtered(end+1)=RR_unfiltered(n);
         T_peak_filtered(end+1)=T_peak(n);
         Count=Count+1;
    end
end

RR = interp1(T_peak_filtered,RR_filtered,TQ,'PCHIP');
RR_unfiltered_ts=interp1(T_peak(1:end-1),RR_unfiltered,TQ,'PCHIP');

if ANOMALIES_DETECTED==1
    figure
    plot(TQ,RR_unfiltered_ts)
    hold;
    plot(TQ,RR,'r--','linewidth',2)
    xlabel('Time (s)')
    ylabel('RR interval')
    legend('RRI with anomalies','Output RRI with anomalies removed')
else
    disp('No anomalies detected.')
    figure
    set(gca, 'FontSize', 16)
    plot(TQ,RR,'r','linewidth',2)
    xlabel('Time (s)')
    ylabel('RR interval')
end    
function [y_peak,T_peak]=ECG_peak_detection_v2(x_ECG,T_ECG,fs,MINPEAKHEIGHT,MINPEAKDISTANCE)
FontSize = 16; 

MINPEAKDISTANCE = ceil(MINPEAKDISTANCE*fs);

% find peaks in ECG
[y_peak,T_peak] = findpeaks(x_ECG,'MINPEAKDISTANCE',MINPEAKDISTANCE,'MINPEAKHEIGHT',MINPEAKHEIGHT);



if isempty(T_peak)
    disp('Warning. No peaks were identified. Adjust the parameters.')
else
    
    number_seconds=ceil(length(x_ECG)/fs);
    num_peaks=length(y_peak);
    bpm=ceil(60*num_peaks/number_seconds);    
    
    if bpm<20
        disp('Warning. A small number of peaks were identified. Recommend adjusting the parameters.')
    end    

    if bpm>300
        disp('Warning. A large number of peaks were identified. Recommend adjusting the parameters.')
    end    

end


T_peak=T_peak/fs;

figure
subplot(2,1,1)
plot(T_ECG,x_ECG,'k')
ylabel('ECG')
hold on;
plot(T_peak,x_ECG(ceil(T_peak*fs)),'markersize',16,'color','r','Marker','.','linestyle','n')
plot(T_ECG,ones(size(x_ECG))*MINPEAKHEIGHT,'b--','linewidth',2)
legend('Filtered ECG','Estimated R peaks','ampthresh')
h1=gca();
subplot(2,1,2)
plot(T_peak(2:end),diff(T_peak),'color','k')
ylabel('Difference between R peaks (s)')
xlabel('Time (s)')
h2=gca();
linkaxes([h1,h2], 'x');


function [xECG,fsECG,hpECG,lpECG,fsRRI,ampthresh,timethresh,anomalyparam,ampthresh_autodetect] = init(varargin)

if nargin<2
    error('At least two inputs must be supplied: ''xECG'' and ''fsECG''.')
end

xECG = varargin{1};
fsECG = varargin{2};


if ~isvector(xECG)
  error('''xECG'' must be a vector.')
end

if size(xECG,1) > 1
  xECG = xECG.';
end

if ~isscalar(fsECG)
  error('''fsECG'' must be a scalar.')
end

if (fsECG)<32
  error('The ECG sampling frequency is insufficiently low for analysis (<32 Hz).')
end


if (fsECG)<120
  disp('Warning. The ECG sampling frequency is below the recommended value (<120 Hz).')
end


if ceil(length(xECG)/fsECG)<10
  error('At least 10 s of ECG are required to calculate RRI.')
end


if nargin>2
    inopts = struct(varargin{3:end});
end
parameter_fields = {'hpECG','lpECG','fsRRI','ampthresh','timethresh','anomalyparam','ampthresh_autodetect'};


% default parameter values
default_params.hpECG = 5;
default_params.lpECG = 20;
default_params.fsRRI = 4;
default_params.ampthresh = 0.5*10^(-3); 
default_params.timethresh = 0.4; 
default_params.anomalyparam = 0.3;
default_params.ampthresh_autodetect = 'Y';
opts = default_params;

if nargin==2
    inopts = default_params;
end

parameters = fieldnames(inopts);

for parameter = parameters'
  if ~any(strcmpi(char(parameter), parameter_fields))
    error(['parameter not recognized: ',char(parameter)])
  end
  if ~isempty(eval(['inopts.',char(parameter)])) 
    eval(['opts.',(char(parameter)),' = inopts.',char(parameter),';'])
  end
end

hpECG = opts.hpECG;
lpECG = opts.lpECG;
fsRRI = opts.fsRRI;
ampthresh = opts.ampthresh;
timethresh = opts.timethresh;
anomalyparam = opts.anomalyparam;
if (any(strcmp(varargin,'ampthresh')))&&(any(strcmp(varargin,'ampthresh_autodetect')))
    disp('Warning. Parameters ''ampthresh'' and ''ampthresh_autodetect'' have BOTH been selected. Auto-detection of the amplitude threshold has been turned OFF.')
    ampthresh_autodetect = 'N';
elseif (any(strcmp(varargin,'ampthresh')))
    ampthresh_autodetect = 'N';
    disp('A parameter value for ''ampthresh'' has been selected. Auto-detection of the amplitude threshold has been turned OFF.')
else
    ampthresh_autodetect = opts.ampthresh_autodetect;
end


if ~ischar(ampthresh_autodetect)
  error('''ampthresh_autodetect'' must be a char.')
elseif (~strcmp(ampthresh_autodetect,'Y'))&&(~strcmp(ampthresh_autodetect,'N'))
    error('''ampthresh_autodetect'' must be ''Y'' or ''N'' ')
end

if ~isscalar(lpECG)
  error('lpECG must be a scalar.')
end

if ~isscalar(hpECG)
  error('hpECG must be a scalar.')
end

if ~isscalar(fsRRI)
  error('fsRRI must be a scalar.')
end

if ~isscalar(ampthresh)
  error('ampthresh must be a scalar.')
end

if ~isscalar(timethresh)
  error('timethresh must be a scalar.')
end

if ~isscalar(anomalyparam)
  error('anomalyparam must be a scalar.')
end
