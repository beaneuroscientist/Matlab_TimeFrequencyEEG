function [pow,itpc,freq,raw]=nwavelet(EEG,varargin)
%% computes wavelet analysis of EEG signal
% Computes wavelets and coherence on a variety of condition and channel
% selections. This does not can store full analysis and computes
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
%                phase angles. This requires a lot of memory and space!!!Default is 'off'.
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
%% Compute Laplacian of data
if strcmpi(inps.Laplace,'on')
    EEG=n_laplacian_nola(EEG); % Warning remove eye channels!
end
%% Create the channels on which to compute the wavelet
[~, labels]=readlocs(EEG.chanlocs);
chidx=[]; % this will contain all of the channels indexes for wavelet
if strcmp(inps.Chans,'all')
    chidx=1:length(labels);
else
    for k=1:length(inps.Chans)
        chidx=union(chidx,find(strcmp(inps.Chans(k),labels)));
    end
end
%% Define the wavelet and compute will be the same across channels
min_freq =  inps.Minfreq;
max_freq = inps.Maxfreq;
num_frex = inps.Nbscales;

% define wavelet parameters
time = -1:1/EEG.srate:1;
if strcmpi(inps.Spacing,'log')
    frex = logspace(log10(min_freq),log10(max_freq),num_frex); %frequencies used logspace
else
    frex = linspace(log10(min_freq),log10(max_freq),num_frex); %frequencies used linear
end
% guasian filter definition
if inps.Cycles==0
    s    = logspace(log10(3),log10(10),num_frex)./(2*pi*frex); %cycles change on a log scale with scale frequency
else
    s = inps.Cycles./(2*pi*frex);
end

%% definte convolution parameters
n_wavelet            = length(time);
n_data               = EEG.pnts*EEG.trials;
n_convolution        = n_wavelet+n_data-1;
n_conv_pow2          = pow2(nextpow2(n_convolution));
half_of_wavelet_size = (n_wavelet-1)/2;

%% Compute the Wavelet iterating over all channels.
%Creates a wavelet variable used for all channel iterations
wavelet=zeros(num_frex,n_conv_pow2);
for fi=1:num_frex
    fprintf('Computing wavelet variable of frequency %d of %d\n',fi, num_frex)
    wavelet(fi,:) = fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 );
end
%% compute baseline index wavelet decomposition on each channel
if isempty(inps.Baseline)
    baseidx=[EEG.times(1),EEG.times(end)];
elseif length(inps.Baseline)==2    
    baseidx = dsearchn(EEG.times',[inps.Baseline(1) inps.Baseline(2)]');
else
    error('Seperate file for baseline is not coded yet')
end

eegfft=zeros(length(chidx),n_conv_pow2);
eegpower = zeros(length(chidx),num_frex,EEG.pnts); % frequencies X time
itpc = zeros(length(chidx),num_frex,EEG.pnts);

if strcmp(inps.RawRetain,'on')
    rawangle=zeros(length(chidx),num_frex,EEG.pnts,EEG.trials);
    raw=zeros(length(chidx),num_frex,EEG.pnts,EEG.trials);
end
if strcmpi(inps.Laplace,'on')   
    data=EEG.laplacian;
else
    data=EEG.data;
end
for ch=1:length(chidx)
    tic;
    fprintf('Computing wavelet on Channel %d\n',ch);
    eegfft(ch,:) = fft(reshape(data(chidx(ch),:,:),1,EEG.pnts*EEG.trials),n_conv_pow2);

    % initialize
    
    for fi=1:num_frex
        eegconv = ifft(wavelet(fi,:).*eegfft(ch,:));
        eegconv = eegconv(1:n_convolution);
        %raw wavelet is below with complex values
        eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
        if strcmp(inps.RawRetain,'on')
            raw(ch,fi,:,:)=reshape(eegconv,EEG.pnts,EEG.trials);
        end
        % compute power
        % temppower = mean(abs(reshape(eegconv,EEG.pnts,EEG.trials)).^2,2); % compute power without conjugation
		% it is computationally faster to compute power using conjugation and it is slightly more accurate implimented below
		% if you take the following code
		% a=rand(20);
		% b=abs(a).^2;
		% c=a.*conj(a);
		% b-c  % you will note that b and c are not indentical after about 15 decimal places for some values.  
		% This is due to an error introduced by the abs function. By taking a square root you are introducing some error
		% alternatively try 
		% temppower = mean(real(reshape(eegconv,EEG.pnts,EEG.trials)).^2+image(reshape(eegconv,EEG.pnts,EEG.trials)).^2,2);
		% Slowest method but more accurate than method abs(A).^2
		temppower = mean(reshape(eegconv,EEG.pnts,EEG.trials).*conj(reshape(eegconv,EEG.pnts,EEG.trials)),2);
        if strcmp(inps.Powmethod,'decibel')
            eegpower(ch,fi,:) = 10*log10(temppower./mean(temppower(baseidx(1):baseidx(2))));
        elseif strcmp(inps.Powmethod,'mean')
            eegpower(ch,fi,:) = temppower./mean(temppower(baseidx(1):baseidx(2)));
        elseif strcmp(inps.Powmethod,'median')
			temppower=median(abs(reshape(eegconv,EEG.pnts,EEG.trials)).^2,2);
            eegpower(ch,fi,:) = 10*log10(temppower./mean(temppower(baseidx(1):baseidx(2))));
        elseif strcmp(inps.Powmethod,'zscore')
             eegpower(ch,fi,:) = (temppowr-mean(temppower(baseidx(1):baseidx(2))))./std(baseidx(1):baseidx(2),0,1);
        end
        % compute Coherenence
        tmp=reshape(eegconv,1,EEG.pnts,EEG.trials);
        tmp=squeeze(tmp);
        tmp=angle(tmp);
        %rawangle(ch,fi,:,:)=tmp;
        for ep=1:EEG.pnts    
            itpc(ch,fi,ep) = abs(mean(exp(1i*tmp(ep,:)),2));
        end
    end
    toc;
end
pow=eegpower;

freq=frex;
%EEG.wavelet.rawangle=rawangle;

