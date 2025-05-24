clear all;
close all;

filename = 'filename.wav';
[x, fs] = audioread(filename);
%sound(x, fs)

% fs = sampling frequency - 44kHz
% x = vector of data
PFA = 10^(-4); 
RefWindow = 32;
guardCells = 2; %Specifies total number of guard cells
RefWindowHalf = RefWindow/2; %Finding number of reference cells on ecah side of CUT
GuardCellsHalf = guardCells/2; %Number of guard cells on each side of CUT
CUT_index = RefWindow/2 + GuardCellsHalf + 1; %Finding index of CUT
sliderLen = RefWindow + guardCells + 1; %Creating a new 'window' to loop signal

%{
Filters for different signals
Fc_Hz = (1.35)*1000; % 5 kHz
Fs = 44.1e3; % Sampling frequency 
w0 = Fc_Hz/(Fs/2);  % Normalised frequency 
bw = 0.005;  % Bandwidth of the notch filter  
[b_notch, a_notch] = iirnotch(w0, bw);
RX_signal_filtered = filtfilt(b_notch, a_notch, x);

Fc_Hz = (0.85)*1000; % 5 kHz
Fs = 44.1e3; % Sampling frequency 
w0 = Fc_Hz/(Fs/2);  % Normalised frequency 
bw = 0.006;  % Bandwidth of the notch filter  
[b_notch, a_notch] = iirnotch(w0, bw);
RX_signal_filtered2 = filtfilt(b_notch, a_notch, RX_signal_filtered);

Fc_Hz = (2.2)*1000; % 5 kHz
Fs = 44.1e3; % Sampling frequency 
w0 = Fc_Hz/(Fs/2);  % Normalised frequency 
bw = 0.003;  % Bandwidth of the notch filter  
[b_notch, a_notch] = iirnotch(w0, bw);
RX_signal_filtered3 = filtfilt(b_notch, a_notch, RX_signal_filtered2);
%}

% Transpose and apply spectrogram on x
y = transpose(x); % Convert x from a column vector to a row vector
FrameLen = 4096; % Value needs to be revised later
ov = 0.5;
[s, f, t] = Sayuri_Spectrogram(y, FrameLen, ov, fs);
[normalised_FFT, v_kmph] =  Sayuri_SubsetNormalisedSpectrogram(s, f);

%{
%Plotting reduced matrix 
figure; imagesc(t, v_kmph, 20*log10(abs(normalised_FFT)));
colormap('jet'); colorbar;                                                                             
xlabel('Time (s)'); 
ylim([0.2 80]);
clim([-60 0]);
ylabel('Speed (km/h)');         
axis xy; 
%}

%Extending the FFT
normalised_FFT_Extended = [normalised_FFT; normalised_FFT([1:sliderLen], :)];

%Power law detector
normalised_FFT_Extended_AfterPowerLawDetector_matrix = abs(normalised_FFT_Extended).^2;

%Finding dimensions of normalised FFT
[NumRows, NumCols] = size(normalised_FFT);

%Finding length of extended matrix
extended_NumRows = NumRows + sliderLen;

totalNoDetections = 0;


%Looping through columns of matrix to apply CA-CFAR
for r = 1:NumCols
    normalised_FFT_column = normalised_FFT_Extended_AfterPowerLawDetector_matrix(:, r);
    [detectCol, numberOfDetections] = Sayuri_CACFAR(RefWindow, guardCells, sliderLen, PFA, normalised_FFT_column, extended_NumRows, CUT_index);
    detectMatrix(:, r) = detectCol;
    totalNoDetections = totalNoDetections + numberOfDetections;
end

% Replacing the first rows that had not been used as CUT with the last rows
% of detectMatrix (so that detect points actually match up) 
detectMatrix([1:(RefWindowHalf+1)], :) = detectMatrix([(NumRows+1):(NumRows+RefWindowHalf+1)], :); 

% Replacing the last rows of extended matrix with zeros because they have
% been copied to beginning of matrix. Otherwise indices higher than NumRows
% are detected which are longer than v_kmph (access error caused)
detectMatrix([NumRows+1:(NumRows+RefWindowHalf+1)], :) = zeros((RefWindowHalf+1), NumCols);
[row, col] = find(detectMatrix == 1);

%Plotting reduced and normalised matrix 
figure; imagesc(t, v_kmph, 20*log10(abs(normalised_FFT)));
colormap('jet'); colorbar;                                                                             
xlabel('Time (s)'); 
ylim([2 60]);
clim([-60 0]);
ylabel('Speed (km/h)');         
axis xy; 


%Plotting reduced and normalised matrix 
t_DetOnly = t(col);
v_kmph_DetOnly = v_kmph(row);

figure; imagesc(t, v_kmph, 20*log10(abs(normalised_FFT)));
colormap('jet'); colorbar;                                                                             
xlabel('Time (s)'); 
ylim([2 60]);
clim([-60 0]);
ylabel('Speed (km/h)');         
axis xy; 

hold on;
p = plot(t_DetOnly, v_kmph_DetOnly, 'kx', 'markersize', 10);
ylim([2 60]);
%p.Color = '#00841a';


%Filtering outliers from detection matrix and finding polynomial model
filtered_v_kmph_DetOnly = filloutliers(v_kmph_DetOnly, "spline");
[q,S] = polyfit(t_DetOnly, filtered_v_kmph_DetOnly, 1);

%Finding length of signal in seconds
N = length(y); % sample length
slength = N/fs;
t_samples = 0:1:round(slength);

%Evaluating polynomial at new time increments
v_poly = polyval(q, t_samples);

%Plotting filtered detections
figure; plot(t_DetOnly, filtered_v_kmph_DetOnly,'x');
xlabel('Time (s)'); 
ylim([2 60]);
ylabel('Speed (km/h)');         
axis xy;

%Plotting new polynomial
figure; plot(t_samples, v_poly, 'o');
legend('Linear Fit');
xlabel('Time (s)');
xlim([0 length(t_samples)]);
ylim([2 60]);
ylabel('Speed (km/h)');
title('Speed vs Time');


function [s, f, t] = Sayuri_Spectrogram(y, FrameLen, ov, fs)
    HopSize = (FrameLen*(1-ov));
    w = transpose(hann(FrameLen));  %Hann window chosen because 
    %Hamming window may a minor discontinuities at edges   
    i = 1; %index to find start of next frame                                    
    FrameEnd = FrameLen; %index of the last frame because frame length 
    % needs to remain constant  as it moves through vector                                                              
    j = 1; %counter                                                                         
    NoFrameOV = ov*FrameLen;                                                           
    while FrameEnd < length(y)                                                          
        f1 = y(i:FrameEnd);                                                       
        i = i + NoFrameOV; %Update start of next frame length                                                                                          
        FrameEnd = FrameEnd + NoFrameOV; %Update end of next frame length                 
        FFT1 = fftshift(fft(f1.*w)); % Apply windowing and FFT                                   
        FFT_matrix(:, j) = transpose(FFT1); %Tranpose and add to matrix                            
        j = j + 1; %Update counter
    end

    f = (-FrameLen/2:1:(FrameLen/2-1))*fs/FrameLen; 
    s = FFT_matrix;
    t = (0:1:(size(FFT_matrix,2)-1))*(1/fs*HopSize);                             
end

function [normalised_FFT, v_kmph_output] = Sayuri_SubsetNormalisedSpectrogram(s, f)
    c = 299792458;
    fc = 24e9;
    lambda = c/fc;

    % Convert f (Hz) to speed (m/s) using Doppler
    v = (f*lambda)/2;

    %Converting (m/s) to (km/h)
    v_kmph = v*3.6;

    %Finding subset of FFT matrix necessary for range of car speeds
    reduced_indices = find(v_kmph <= 70 & v_kmph >= 1.6);
    reduced_FFT = s(reduced_indices, :);

    v_kmph_output = v_kmph(reduced_indices);

    %Finding the maximum value of FFT matrix
    max_FFT_Value = max(max((reduced_FFT)));

    %Normalising and plotting reduced matrix 
    normalised_FFT = reduced_FFT/max_FFT_Value;     
end

function [detectCol, numberOfDetections] = Sayuri_CACFAR(RefWindow, guardCells, sliderLen, PFA, DataAfterPowerLawDetector_column, extended_NumRows, CUT_index)
    RefWindowHalf = RefWindow/2; %Finding number of reference cells on each side of CUT
    i = 1; %index to find start of next frame
    indexBeforeExcludedCells = RefWindowHalf; %To find index before
    %guard cells and CUT
    indexAfterExcludedCells = RefWindowHalf + guardCells + 1 + 1; % To find index
    %after guard cells and CUT
    sliderEnd = sliderLen; %Last index of slider

    alphaCA = ((PFA^(-1/RefWindow)) - 1); %CFAR constant
    gCA = [];
    TCA = []; %CA-CFAR threshold

   while sliderEnd < extended_NumRows   
        gCA(i, :) = sum(DataAfterPowerLawDetector_column(i:indexBeforeExcludedCells)) + sum(DataAfterPowerLawDetector_column(indexAfterExcludedCells:sliderEnd));   
        TCA(i, :) = gCA(i, :)*alphaCA; 
        i = i + 1; %Update start of next frame length    
        indexBeforeExcludedCells = indexBeforeExcludedCells + 1;       
        indexAfterExcludedCells = indexAfterExcludedCells + 1;      
        sliderEnd = sliderEnd + 1; %Update end of next frame length     
   end
   NumRows = extended_NumRows - sliderLen;
   numberOfDetections = 0;
   counter = 1;
   detectCol = [];
   % Looping through the real signal and comparing it to the detecion threshold
   % Comparing values from CUT (index 18) to TCA (index 1)
   for j = CUT_index:(NumRows + RefWindowHalf + 1)
       % Adding 1 or 0 based on whether there is a detection or not
       if (DataAfterPowerLawDetector_column(j) >= TCA(counter)) 
            detectCol(j) = 1;
            numberOfDetections = numberOfDetections + 1;
            counter = counter + 1;
       else
            detectCol(j) = 0;
            counter = counter + 1;
       end
   end
   % Storing indices where detections are found
   detect_indices = find(detectCol == 1);
end