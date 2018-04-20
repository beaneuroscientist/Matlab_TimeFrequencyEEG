function EEG=nwavelet_padding(EEG,varargin)
%% computes wavelet analysis of EEG signal by padding the data set
% Computes wavelets and coherence on a variety of condition and channel
% selections. This does not store full analysis and computes
% averages across conditions 
%  EEG-  EEG data structure
%  Optional inputs
%  'Chans'- cell array of channel names to compute wavelet on.  Default
%           'all' computes on all channels present.
%  'Cycles'- number of wavelet cycles or if set to 0 will change cycles
%            with frequencies from 3-10 cylces. Will generate a warning if cylces is less 
%            than 3 or greater than 14. Default 6. 
%  'Laplace'- 'on' or 'off' use surface laplacian instead of EEG.data.
%              Will compute using Perin method if laplacian is not computed.
%              Default 'on'
%  'Minfreq'- min frequency to use default is 3
%  'Maxfreq'- max frequency to use. default is 100
%  'Nbscales'- number of frequency scales to use.  Default is 100.
%  'Spacing'- 'linear' or 'log' spacing of scales as log or linear spaced
%              default is log spacing
%  'Baseline'- array.  if [] then will use all points and trials. If of
%                      form with two numbers [-200 0] will use time points
%                      from EEG as baseline, for example use average
%                      power on baseline.  If the array is longer than 2
%                      then will compute wavelet power on array provided
%                      and use as the baseline. default is all times before
%                      zero [min(EEG.times) 0]
%  'Powmethod'- method for computing baseline of power.  'mean' subtracts
%               mean from baseline from all power points. 'median'
%               subtracts median of data from all power points.  'decibel'
%               uses a decibel computation for power.  'zscore' uses a
%               zscore for power computation. 
%  'Condgroups'- matrix or cell array of condition groupings for power and
%                coherence analysis where each row is a grouping of
%                conditions. for example [1 2; 3 4] will create two
%                seperate powers based on the groupings of epochs 1 and 2;
%                and 3 and 4.  Similar for cell arrays when epochs are
%                strings. Default is average all epochs.
%  'Condlabels'- required if condgroups is set to none default. Label names
%                for each row of condlabels that is set. For example
%                {'rare','frequent'}. EEG.wavelet.power.rare and *.frequent
%                will be created with average across all trials in each group.
%  'RawRetain'- 'on' or 'off' retain raw wavelet across all trials and
%                phase angles. This requires a lot of memory and space!!!Default is 'on'.
%% create the inputs

try
inps=finputcheck(varargin, ...
               { 'Chans'    'cell' [] {'all'};
                 'Cycles'   'integer'   [0 1000]       6;
                 'Laplace'  'string'   {'on', 'off'}, 'off';       
                 'Minfreq'  'real'     [0.0001 1000]   3;
                 'Maxfreq'  'real'     [1 10000]       100;
                 'Nbscales' 'integer'  [1 10000]       100;
                 'Spacing'  'string'   {'linear', 'log'}, 'log';         
                 'Baseline' 'real'     []           [min(EEG.times) 0];
                 'Powmethod' 'string'  {'mean', 'median', 'decibel','zscore'} 'decibel';
                 'Condgroups' 'real'   []   [];
                 'Condlabels' 'cell'   []   {'all'};
                 'RawRetain'  'string'   {'on', 'off'}, 'off';});
catch
    inps=finputcheck(varargin, ...
               { 'Chans'    'cell' [] {'all'};
                 'Cycles'   'integer'   [0 1000]       6;
                 'Laplace'  'string'   {'on', 'off'}, 'off';       
                 'Minfreq'  'real'     [0.0001 1000]   3;
                 'Maxfreq'  'real'     [1 10000]       100;
                 'Nbscales' 'integer'  [1 10000]       100;
                 'Spacing'  'string'   {'linear', 'log'}, 'log';         
                 'Baseline' 'real'     []           [min(EEG.times) 0];
                 'Powmethod' 'string'  {'mean', 'median', 'decibel','zscore'} 'decibel';
                 'Condgroups' 'cell'   []   [];
                 'Condlabels' 'cell'   []   {'all'};
                 'RawRetain'  'string'   {'on', 'off'}, 'off';});
end
times=EEG.times;
len=length(EEG.times);
EEG.data=padarray(EEG.data,[0,2*len,0],'symmetric','both');
delT=abs(diff(times));
t1=min(EEG.times);
t2=max(EEG.times);
times2=t1-4*2*len:delT:t2+4*2*len;
EEG.pnts=5*len;
EEG.times=times2;
if strcmp(inps.RawRetain,'on')
[pow,itpc,freq,raw]=nwavelet(EEG, inps);
EEG.wavelet.raw=raw(:,:,2*len+1:3*len,:);
else
    [pow,itpc,freq]=nwavelet(EEG, inps);
end

%% unpad the data
EEG.data=EEG.data(:,2*len+1:3*len,:);
EEG.wavelet.power=pow(:,:,2*len+1:3*len);
EEG.wavelet.itpc=itpc(:,:,2*len+1:3*len);

EEG.times=t1:delT:t2;
EEG.pnts=len;
EEG.frequencies=freq;
