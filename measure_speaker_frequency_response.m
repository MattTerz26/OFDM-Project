% measure_speaker_frequency_response.m
% Spielt einen logarithmischen Sweep ab, zeichnet das Mikrofonsignal auf
% und schätzt daraus die Frequenzgang-Magnitude und Phase des Lautsprechers.

clear; close all; clc;

%% Benutzer-Parameter
fs = 48e3;                 % Abtastrate [Hz]
sweepDuration = 5;         % Dauer des Sweeps [s]
startSilence = 0.5;        % Stille vor dem Sweep [s]
endSilence = 0.5;          % Stille nach dem Sweep [s]
playbackLevel = 0.6;       % 0..1, Lautstärke für den Sweep
micDeviceID = [];          % ID des Aufnahmegeräts (leer = Standard)
speakerDeviceID = [];      % ID des Wiedergabegeräts (leer = Standard)

%% Sweep generieren (logarithmisch von 20 Hz bis Nyquist)
t = (0:1/fs:sweepDuration).';
sweep = chirp(t, 20, sweepDuration, fs/2, 'logarithmic');
sweep = playbackLevel * sweep;
signalToPlay = [zeros(round(startSilence * fs), 1); sweep; zeros(round(endSilence * fs), 1)];
totalDuration = numel(signalToPlay) / fs;

%% Aufnahme vorbereiten
bitsPerSample = 24;
if isempty(micDeviceID)
    recorder = audiorecorder(fs, bitsPerSample, 1); % Standard-Mikro
else
    recorder = audiorecorder(fs, bitsPerSample, 1, micDeviceID);
end

if isempty(speakerDeviceID)
    player = audioplayer(signalToPlay, fs, bitsPerSample); % Standard-Ausgabe
else
    player = audioplayer(signalToPlay, fs, bitsPerSample, speakerDeviceID);
end

%% Abspielen und gleichzeitig aufnehmen
record(recorder);          % Aufnahme starten
pause(0.1);                % kleiner Vorlauf für stabilen Start
play(player);              % Sweep abspielen
pause(totalDuration + 0.5);% warten bis alles abgespielt wurde
stop(recorder);

%% Aufgenommene Daten auslesen und auf Länge des Sweeps kürzen
micSignal = getaudiodata(recorder);
minLen = min(numel(signalToPlay), numel(micSignal));
signalToPlay = signalToPlay(1:minLen);
micSignal = micSignal(1:minLen);

%% Frequenzgang schätzen (einfache FFT-Division)
nfft = 2^nextpow2(minLen);
X = fft(signalToPlay, nfft);
Y = fft(micSignal, nfft);
H = Y ./ (X + eps);                      % eps verhindert Division durch 0
freqs = (0:nfft/2).' * (fs / nfft);
H_half = H(1:numel(freqs));
magnitudeDb = 20 * log10(abs(H_half));
phaseDeg = unwrap(angle(H_half)) * 180/pi;

%% Linearsten 2-kHz-Bereich in der Magnitude finden
targetBW = 2000; % Hz
maxFreq = fs/2;
bestErr = inf;
bestRange = [NaN NaN];
bestMask = false(size(freqs));

for idx = 1:numel(freqs)
    fStart = freqs(idx);
    fEnd = fStart + targetBW;
    if fEnd > maxFreq
        break; % restliche Startpunkte wären ebenfalls zu hoch
    end
    mask = freqs >= fStart & freqs <= fEnd;
    if nnz(mask) < 4
        continue; % zu wenige Punkte für Regression
    end
    p = polyfit(freqs(mask), magnitudeDb(mask), 1); % lineare Regression
    fitVals = polyval(p, freqs(mask));
    err = sqrt(mean((fitVals - magnitudeDb(mask)).^2)); % RMS-Fehler
    if err < bestErr
        bestErr = err;
        bestRange = [fStart, fEnd];
        bestMask = mask;
    end
end

bestCenterFreq = mean(bestRange);

%% Plot erstellen
figure('Name', 'Speaker Frequency Response', 'Position', [200 200 900 600]);
subplot(2,1,1);
semilogx(freqs, magnitudeDb, 'LineWidth', 1.3);
hold on;
if ~isnan(bestCenterFreq)
    ylims = ylim;
    fill([bestRange(1) bestRange(2) bestRange(2) bestRange(1)], ...
         [ylims(1) ylims(1) ylims(2) ylims(2)], ...
         [0.85 0.9 1], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    semilogx(freqs(bestMask), magnitudeDb(bestMask), 'r', 'LineWidth', 1.6);
    ylim(ylims); % Rücksetzen, falls fill geändert hat
end
grid on; xlim([20 fs/2]);
xlabel('Frequency [Hz]'); ylabel('Magnitude [dB]');
title('Estimated Frequency Response');

subplot(2,1,2);
semilogx(freqs, phaseDeg, 'LineWidth', 1.0);
grid on; xlim([20 fs/2]);
xlabel('Frequency [Hz]'); ylabel('Phase [deg]');

%% Plot speichern
plotDir = fullfile(pwd, 'plots_frequency_response');
if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
plotPath = fullfile(plotDir, ['frequency_response_' timestamp '.png']);
saveas(gcf, plotPath);

fprintf('Frequenzgang geplottet und gespeichert unter: %s\n', plotPath);
fprintf('Hinweis: Falls andere Geräte genutzt werden sollen, "audiodevinfo" ausführen und IDs oben setzen.\n');
if ~isnan(bestCenterFreq)
    fprintf('Linearster 2-kHz-Bereich: %.1f Hz - %.1f Hz, Center: %.1f Hz, RMS-Fehler: %.2f dB\n', ...
        bestRange(1), bestRange(2), bestCenterFreq, bestErr);
else
    fprintf('Kein gültiger Bereich für die Linearisierung gefunden.\n');
end
