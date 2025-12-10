%% Diagonal Checkerboard Pilot/Data Matrix Plot (Standalone)

clear; clc;

%% Manual parameters
pilotSpacing = 140;   % pilot spacing in frequency
pilotShift   = 1;     % shift per OFDM symbol
N            = 512;   % number of subcarriers
nDataOfdmSym = 50;    % number of data symbols (1st = training)

%% Pilot index computation
baseIdx = 1:pilotSpacing:N;

pilotMap = zeros(N, nDataOfdmSym + 1);

pilotMap(:,1) = 2; % full training symbol

for s = 1:nDataOfdmSym
    col   = s + 1;
    shift = (s-1) * pilotShift;

    pilotIdx_s = mod(baseIdx - 1 + shift, N) + 1;
    pilotIdx_s = sort(pilotIdx_s);

    pilotMap(pilotIdx_s, col) = 1;
end

%% Plot
figure;
imagesc(pilotMap);
colormap([
    0 0 0;        % data
    0 0.7 1;      % pilot
    1 0.3 0.3     % training
]);
colorbar('Ticks',[0,1,2], 'TickLabels',{'Data','Pilot','Training'});

xlabel('OFDM Symbol Index');
ylabel('Subcarrier Index');
title(sprintf('Diagonal Pilot/Data Grid (spacing = %d, shift = %d)', ...
              pilotSpacing, pilotShift));

set(gca,'YDir','normal');

%% Save plot in same folder as script
scriptPath = mfilename('fullpath');
scriptDir  = fileparts(scriptPath);

filename = sprintf('pilotGrid_spacing%d_shift%d.png', pilotSpacing, pilotShift);
savePath = fullfile(scriptDir, filename);

saveas(gcf, savePath);

fprintf('Saved plot to: %s\n', savePath);
