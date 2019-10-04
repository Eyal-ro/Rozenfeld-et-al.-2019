function [v Isyn spike_times spikeTrain spikes_per_sec] = LIF_model(fr,g1)
%% PN model function

%   Adapted from: "Nagel, K. I., Hong, E. J., & Wilson, R. I. (2015).
%   Synaptic and circuit mechanisms promoting broadband transmission of olfactory stimulus dynamics. Nature neuroscience"
%
% [v Isyn spike_times spikeTrain spikes_per_sec] = LIF_PNmodel_many_inputs(fr,g1)
%   Input:
%       fr: A column vector of ORN firing rates (in spikes/sec)
%       g1: conductance gain (keep at 1)
%   Output:
%       v: membrane voltage
%       Isyn: synaptic current
%       spike_times: spike times (sampels unit)
%       spikeTrain: binary spike train vector
%       spikes_per_sec: PSTH

s = fr/1000;        % spike rate in spikes/ms



%% depression model parameters
f1 = 0.23;       % depletion rate for fast component
tauD1 = 1006;   % depletion recovery time (ms)
% g1 = 20;      % conductance gain
tauC1 = 9.3;   % conductance decay time constant (ms)

f2 = 1  %0.0073;       % depletion rate for slow component
tauD2 = 33247;   % depletion recovery time (ms)
g2 = 1.8;      % conductance gain
tauC2 = 80;   % conductance decay time constant (ms)

%% depression model (run at 1KHz, time step is 1ms)
for j=1:size(s,2)
A1(1,j) = 1./(tauD1*f1*s(1,j)+1);   % initial A1
c1(1,j) = A1(1,j)*s(1,j)*g1*tauC1;    % initial c1

A2(1,j) = 1./(tauD2*f2*s(1,j)+1);   % initial A2
c2(1,j) = A2(1,j)*s(1,j)*g2*tauC2;    % initial c2

for i=1:length(s)-1
    
    A1(i+1,j) = A1(i,j) + (-f1*s(i,j)*A1(i,j) + (1-A1(i,j))/tauD1); % depletion
    c1(i+1,j) = c1(i,j) + (g1*A1(i,j)*s(i,j) - c1(i,j)/tauC1);  % conductance
    
    A2(i+1,j) = A2(i,j) + (-f2*s(i,j)*A2(i,j) + (1-A2(i,j))/tauD2); % depletion
    c2(i+1,j) = c2(i,j) + (g2*A2(i,j)*s(i,j) - c2(i,j)/tauC2);  % conductance
    
end
end

%% membrane potential model parameters
rm = 800/1000;             % input resistance (GOhms)
Esyn = -10;                % excitatory synaptic reversal potential (mv)
Einh = -70;                % inhibitory reversal potential (mv)
Eleak = -70;               % leak reversal potential (mv)
tauM = 5;                  % membrane time constant (ms)
thresh = -40;
h = 0.1;                    % time step (this part of the model is run at 10kHz)
t(1) = 0;                   % time vector
v(1) = -70;                 % initial membrane voltage
arp = 2;                    % this being 0.02 doesn't work.
relative_refrac = -10;      % mV

for j=1:size(s,2)
gsyn(:,j) = resample(c1(:,j)+c2(:,j),10,1);            % excitatory synaptic conductance in nS; upsample to 10kHz
% gsyn(:,j) = gsyn(:,j)*0.7; % scale synaptic conductance by constant

ginh = zeros(size(gsyn));               % inhibitory synaptic conductance (can be used for post-synaptic inhibition)
end

gsyn=sum(gsyn,2);
% gsyn = gsyn * 0.7; % scale synaptic conductance by constant
spike_times(1)=0;
num_spikes = 0;

% membrane potential model
for i=1:length(gsyn)
    
    Isyn(i) = gsyn(i)*rm*(v(i)-Esyn);       %  excitatory synaptic current
    v(i+1) = v(i) - h/tauM*(v(i)-Eleak + gsyn(i)*rm*(v(i)-Esyn) + ginh(i)*rm*(v(i)-Einh));      %  membrane potential
    t(i+1) = t(i) + h;                      % time vector
    %----------------------
    % spike
    if (v(i+1) > thresh) 
        if num_spikes>0
            % check we're not in absolute refac period
            if (t(i)>=(spike_times(end)+arp))
                % reset voltage
                v(i+1) = Eleak + relative_refrac;

                v(i) = -0;
                % record spike
                spike_times(end+1) = t(i);
                % increment spike count
                num_spikes = num_spikes+1;
            end
        else
            % no spikes yet - no need to check
            v(i+1) = Eleak;

            v(i) = -0;
            % record spike
            spike_times(end+1) = t(i);
            % increment spike count
            num_spikes = num_spikes+1;
        end
    end
    %----------------------

    
    
end

v(end) = [];       % trim vectors to same length
t(end) = [];

% plot results
spikeTrain = zeros(1,round(t(end)));
spikeTrain(round(spike_times(2:end))) = 1; %
% PSTH
binsize = 1000;
lastBin = round(t(end));
edges = 0 : binsize : lastBin;	
x = (mod(spike_times-1,lastBin)+1);
r = (histc(x,edges)*1000) / (binsize);
spikes_per_sec = r(1:end-1);

