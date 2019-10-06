%% Run each section to generate the designated figure pannel
% The code was written with MATLAB R2018a and may not be compatible with previous versions

%% Figure 7B

% For the entire script:
%   fr: input firing rate (spikes/sec)
%   Isyn: synaptic current
%   spike_times: the time of action potentials.
%   Isyn: synaptic current

clc;clear;
PNs = struct('fr', [],'Isyn', [],'spike_times', [],'spikes_per_sec', [],'spikeTrain', [],'Vm', []); % Structure contacting the data
LNs = struct('fr', [],'Isyn', [],'spike_times', [],'spikes_per_sec', [],'spikeTrain', [],'Vm', []); % Structure contacting the data


for k=1:100 % Number of cells
    
    clc;clearvars -except k PNs LNs;
    number_of_ORNs = 20;
    jitter = 10; % ms
    baseline_firing_rate = 8; % spikes/sec
    PNs_ORNs_firing_rates = 140; % spikes/sec
    LNs_ORNs_firing_rates = 140:-5:45; % spikes/sec
    g1 = 1;      % conductance gain
    
    for i=LNs_ORNs_firing_rates
        
        temp_index = double(LNs_ORNs_firing_rates==i);
        a = find(temp_index==1);
        demo_pulse = zeros(500+ randi(jitter),1)+baseline_firing_rate; demo_pulse(200:225)  = i;
        fr_LNs_pulse = repmat(demo_pulse,[40,1]);
        fr_LNs(:,a) = zeros(60000,1)+baseline_firing_rate;
        fr_LNs(20000 : length(fr_LNs_pulse)+19999,a) = fr_LNs_pulse;
        
    end
    
    [v_LNs Isyn_LNs spike_times_LNs spikeTrain_LNs spikes_per_sec_LNs] = LIF_model(fr_LNs,g1);
    
    
    clc;clearvars -except v_LNs Isyn_LNs spike_times_LNs spikeTrain_LNs spikes_per_sec_LNs fr_LNs k PNs LNs LNs_ORNs_firing_rates;
    number_of_ORNs = 20;
    jitter = 10; % ms
    baseline_firing_rate = 8; % spikes/sec
    PNs_ORNs_firing_rates = 140; % spikes/sec
    LNs_ORNs_firing_rates = 140:-5:45; % spikes/sec
    g1 = 1;      % conductance gain
    
    
    for j=1:number_of_ORNs
        
        demo_pulse = zeros(500+ randi(jitter),1)+baseline_firing_rate; demo_pulse(200:225)  = PNs_ORNs_firing_rates;
        fr_PNs_pulse = repmat(demo_pulse,[40,1]);
        fr_PNs(:,j) = zeros(60000,1)+baseline_firing_rate;
        fr_PNs(20000 : length(fr_PNs_pulse)+19999,j) = fr_PNs_pulse;
        
    end
    [v_PNs Isyn_PNs spike_times_PNs spikeTrain_PNs spikes_per_sec_PNs] = LIF_model(fr_PNs,g1);
    
    figure;
    a = area(spikes_per_sec_PNs);a.FaceAlpha = 0.5;
    hold on; b = area(spikes_per_sec_LNs); b.FaceAlpha = 0.5;
    title('Raw PSTH');ylabel('Spikes/sec');
    
    figure;
    a = area(spikes_per_sec_PNs/max(spikes_per_sec_PNs));a.FaceAlpha = 0.5;
    hold on; b = area(spikes_per_sec_LNs/max(spikes_per_sec_LNs)); b.FaceAlpha = 0.5;
    title('Normalized PSTH');ylabel('Spikes/sec');
    
    PNs(k).fr = fr_PNs; PNs(k).Isyn = Isyn_PNs; PNs(k).spike_times = spike_times_PNs;
    PNs(k).spikes_per_sec = spikes_per_sec_PNs; PNs(k).spikeTrain = spikeTrain_PNs; PNs(k).Vm = v_PNs;
    
    LNs(k).fr = fr_LNs; LNs(k).Isyn = Isyn_LNs; LNs(k).spike_times = spike_times_LNs;
    LNs(k).spikes_per_sec = spikes_per_sec_LNs; LNs(k).spikeTrain = spikeTrain_LNs; LNs(k).Vm = v_LNs;
    
    disp(k);
    close all;
end

clc;clearvars -except PNs LNs number_of_ORNs PNs_ORNs_firing_rates LNs_ORNs_firing_rates;


for j=1:length(PNs)
    PNs_PSTH_mat(j,:) = PNs(j).spikes_per_sec;
end

for j=1:length(LNs)
    LNs_PSTH_mat(j,:) = LNs(j).spikes_per_sec;
end

for j=1:length(PNs)
    PNs_current_mat(j,:) = PNs(j).Isyn;
end

for j=1:length(LNs)
    LNs_current_mat(j,:) = LNs(j).Isyn;
end


figure;
subplot(3,1,1);a = area(mean(PNs_PSTH_mat)/max(mean(PNs_PSTH_mat)));a.FaceAlpha = 0.5;box off;set(gca,'TickDir','out');
hold on; subplot(3,1,1);b = area(mean(LNs_PSTH_mat)/max(mean(LNs_PSTH_mat))); b.FaceAlpha = 0.5;box off;set(gca,'TickDir','out');
title('Normalized PSTH');ylabel('Spikes/sec (Norm)');xlabel('Time (sec)');
suptitle({'Simulation results - 100 Runs',strcat('Number of ORNs: ',num2str(number_of_ORNs)),...
    strcat('ORNs firing rate for PNs:',num2str(PNs_ORNs_firing_rates)),'ORNs firing rate for LNs: 140:-5:45'});




for j=1:length(PNs)
    PNs_spikes_mat(j,:) = PNs(j).spikeTrain;
end

for j=1:length(LNs)
    LNs_spikes_mat(j,:) = LNs(j).spikeTrain;
end

subplot(3,1,2);plotSpikeRaster(logical(PNs_spikes_mat),'PlotType','scatter');ylabel('Cells');xlabel('Time (sec)');title('PNs raster plot');box off;set(gca,'TickDir','out');

subplot(3,1,3);plotSpikeRaster(logical(LNs_spikes_mat),'PlotType','scatter');ylabel('Cells');xlabel('Time (sec)');title('LNs raster plot');box off;set(gca,'TickDir','out');

%% Figure 7C

clc;clear;
PNs = struct('fr', [],'Isyn', [],'spike_times', [],'spikes_per_sec', [],'spikeTrain', [],'Vm', []);
LNs = struct('fr', [],'Isyn', [],'spike_times', [],'spikes_per_sec', [],'spikeTrain', [],'Vm', []);


for k=1:100
    
    clc;clearvars -except k PNs LNs;
    number_of_ORNs = 20;
    jitter = 10; % ms
    baseline_firing_rate = 8; % spikes/sec
    PNs_ORNs_firing_rates = 140; % spikes/sec
    LNs_ORNs_firing_rates = 180:-5:85; % spikes/sec
    g1 = 1;      % conductance gain
    
    for i=LNs_ORNs_firing_rates
        temp_index = double(LNs_ORNs_firing_rates==i);
        a = find(temp_index==1);
        demo_pulse = zeros(500+ randi(jitter),1)+baseline_firing_rate; demo_pulse(200:225)  = i;
        fr_LNs_pulse = repmat(demo_pulse,[40,1]);
        fr_LNs(:,a) = zeros(60000,1)+baseline_firing_rate;
        fr_LNs(20000 : length(fr_LNs_pulse)+19999,a) = fr_LNs_pulse;
        
    end
    
    [v_LNs Isyn_LNs spike_times_LNs spikeTrain_LNs spikes_per_sec_LNs] = LIF_model(fr_LNs,g1);
    
    
    clc;clearvars -except v_LNs Isyn_LNs spike_times_LNs spikeTrain_LNs spikes_per_sec_LNs fr_LNs k PNs LNs LNs_ORNs_firing_rates;
    number_of_ORNs = 20;
    jitter = 10; % ms
    baseline_firing_rate = 8; % spikes/sec
    PNs_ORNs_firing_rates = 140; % spikes/sec
    LNs_ORNs_firing_rates = 180:-5:85; % spikes/sec
    g1 = 1;      % conductance gain
    
    
    for j=1:number_of_ORNs
        
        demo_pulse = zeros(500+ randi(jitter),1)+baseline_firing_rate; demo_pulse(200:225)  = PNs_ORNs_firing_rates;
        fr_PNs_pulse = repmat(demo_pulse,[40,1]);
        fr_PNs(:,j) = zeros(60000,1)+baseline_firing_rate;
        fr_PNs(20000 : length(fr_PNs_pulse)+19999,j) = fr_PNs_pulse;
        
    end
    [v_PNs Isyn_PNs spike_times_PNs spikeTrain_PNs spikes_per_sec_PNs] = LIF_model(fr_PNs,g1);
    
    figure;
    a = area(spikes_per_sec_PNs);a.FaceAlpha = 0.5;
    hold on; b = area(spikes_per_sec_LNs); b.FaceAlpha = 0.5;
    title('Raw PSTH');ylabel('Spikes/sec');
    
    figure;
    a = area(spikes_per_sec_PNs/max(spikes_per_sec_PNs));a.FaceAlpha = 0.5;
    hold on; b = area(spikes_per_sec_LNs/max(spikes_per_sec_LNs)); b.FaceAlpha = 0.5;
    title('Normalized PSTH');ylabel('Spikes/sec');
    
    PNs(k).fr = fr_PNs; PNs(k).Isyn = Isyn_PNs; PNs(k).spike_times = spike_times_PNs;
    PNs(k).spikes_per_sec = spikes_per_sec_PNs; PNs(k).spikeTrain = spikeTrain_PNs; PNs(k).Vm = v_PNs;
    
    LNs(k).fr = fr_LNs; LNs(k).Isyn = Isyn_LNs; LNs(k).spike_times = spike_times_LNs;
    LNs(k).spikes_per_sec = spikes_per_sec_LNs; LNs(k).spikeTrain = spikeTrain_LNs; LNs(k).Vm = v_LNs;
    
    disp(k);
    close all;
end

clc;clearvars -except PNs LNs number_of_ORNs PNs_ORNs_firing_rates LNs_ORNs_firing_rates;


for j=1:length(PNs)
    PNs_PSTH_mat(j,:) = PNs(j).spikes_per_sec;
end

for j=1:length(LNs)
    LNs_PSTH_mat(j,:) = LNs(j).spikes_per_sec;
end

for j=1:length(PNs)
    PNs_current_mat(j,:) = PNs(j).Isyn;
end

for j=1:length(LNs)
    LNs_current_mat(j,:) = LNs(j).Isyn;
end


figure;
subplot(3,1,1);a = area(mean(PNs_PSTH_mat)/max(mean(PNs_PSTH_mat)));a.FaceAlpha = 0.5;box off;set(gca,'TickDir','out');
hold on; subplot(3,1,1);b = area(mean(LNs_PSTH_mat)/max(mean(LNs_PSTH_mat))); b.FaceAlpha = 0.5;box off;set(gca,'TickDir','out');
title('Normalized PSTH');ylabel('Spikes/sec (Norm)');xlabel('Time (sec)');
suptitle({'Simulation results - 100 Runs',strcat('Number of ORNs: ',num2str(number_of_ORNs)),...
    strcat('ORNs firing rate for PNs:',num2str(PNs_ORNs_firing_rates)),'ORNs firing rate for LNs: 180:-5:85'});




for j=1:length(PNs)
    PNs_spikes_mat(j,:) = PNs(j).spikeTrain;
end

for j=1:length(LNs)
    LNs_spikes_mat(j,:) = LNs(j).spikeTrain;
end

subplot(3,1,2);plotSpikeRaster(logical(PNs_spikes_mat),'PlotType','scatter');ylabel('Cells');xlabel('Time (sec)');title('PNs raster plot');box off;set(gca,'TickDir','out');

subplot(3,1,3);plotSpikeRaster(logical(LNs_spikes_mat),'PlotType','scatter');ylabel('Cells');xlabel('Time (sec)');title('LNs raster plot');box off;set(gca,'TickDir','out');

%% Figure 7D
clc;clear; close all

LNs = struct('spikes_per_sec', [],'jitter',[]);


for k=1:100
    
    clc;clearvars -except k PNs LNs;
    number_of_ORNs = 20;
    jitter = 10:5:30; % ms
    baseline_firing_rate = 8; % spikes/sec
    LNs_ORNs_firing_rates = 140%180:-5:85; % spikes/sec
    g1 = 1;      % conductance gain
    
    for jit = jitter
        for j=1:number_of_ORNs
            
            demo_pulse = zeros(500+ randi(jit),1)+baseline_firing_rate; demo_pulse(200:225)  = LNs_ORNs_firing_rates;
            fr_LNs_pulse = repmat(demo_pulse,[40,1]);
            fr_LNs(:,j) = zeros(60000,1)+baseline_firing_rate;
            fr_LNs(20000 : length(fr_LNs_pulse)+19999,j) = fr_LNs_pulse;
            
            
        end
        
        [v_LNs Isyn_LNs spike_times_LNs spikeTrain_LNs spikes_per_sec_LNs] = LIF_model(fr_LNs,g1);
        
        
        LNs(end+1) = struct('spikes_per_sec', [spikes_per_sec_LNs],'jitter',[jit])
        
        
    end
    
    disp(k);
    close all;
end

LNs(1) = [];

for j=1:length(LNs)
    LNs_PSTH_mat(j,:) = LNs(j).spikes_per_sec;
    LNs_PSTH_mat(j,:) = LNs_PSTH_mat(j,:)/max(LNs_PSTH_mat(j,:));
end


Jit_10 = mean(LNs_PSTH_mat(1:9:end,:));
Jit_30 = mean(LNs_PSTH_mat(5:5:end,:));

figure;
a = area(Jit_10);a.FaceAlpha = 0.5;box off;set(gca,'TickDir','out');
hold on;
a = area(Jit_30);a.FaceAlpha = 0.5;box off;set(gca,'TickDir','out');

title('Normalized PSTH');ylabel('Spikes/sec (Norm)');xlabel('Time (sec)');

