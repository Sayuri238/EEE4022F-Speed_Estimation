clear all;
close all;

% Constant
fc = 24e9;        % IPM365 has a centre frequency of 24 GHz
c = 299792458;    % speed of light 
R0 = 100;         % Initial distance of vehicle. in meters
v_constant = 11;   % Speed of the vehicle. In meter/second
fs = 44e3;        % Sampling frequency of ADC. 44 kHz
SimulationTime_s = 30;  % Measurement duration
H = 1.5; %Height of radar

PFA = 10^(-3); 
RefWindow = 32;
guardCells = 2; %Specifies total number of guard cells
RefWindowHalf = RefWindow/2; %Finding number of reference cells on ecah side of CUT
GuardCellsHalf = guardCells/2; %Number of guard cells on each side of CUT
CUT_index = RefWindow/2 + GuardCellsHalf + 1; %Finding index of CUT
sliderLen = RefWindow + guardCells + 1; %Creating a new 'window' to loop signal

% Computations
M = SimulationTime_s*fs;  % M is total number of samples in our measurement
lambda = c/fc;
T = (0:1:(M-1))*1/fs;

% Generate received signal
R = sqrt((R0 - v_constant*T).^2 + H^2); %Adjusted R(t) equation to account for change in height
Phi = 4*pi/lambda*(R);     % Phase 
y = sin(Phi);                    % Real values captured from IPM365, from soundcard 

% Plot target range versus time
figure; 
plot(T, R)
xlabel("Time (s)");
ylabel("Distance (m)")
grid on;

% Transpose and apply spectrogram on x
FrameLen = 4096;
ov = 0.5;
[s, f, t] = Sayuri_Spectrogram(y, FrameLen, ov, fs);

%Finding normalised FFT
[normalised_FFT, v_kmph] =  Sayuri_SubsetNormalisedSpectrogram(s, f);

%Extending matrix by sliderLen to wrap around signal
normalised_FFT_Extended = [normalised_FFT; normalised_FFT([1:sliderLen], :)];

%Power law detector
normalised_FFT_Extended_AfterPowerLawDetector_matrix = abs(normalised_FFT_Extended).^2;

%Finding dimensions of normalised FFT
[NumRows, NumCols] = size(normalised_FFT);

%Finding length of extended matrix
extended_NumRows = NumRows + sliderLen;

%Plotting reduced and normalised matrix  
figure; imagesc(t, v_kmph, 20*log10(abs(normalised_FFT)));
colormap('jet'); colorbar;                                                                             
xlabel('Time (s)'); 
ylim([2 50]);
clim([-60 0]);
ylabel('Speed (km/h)');         
axis xy;

detectMatrix = [];
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

%Plotting reduced and normalised matrix and detections
t_DetOnly = t(col);
v_kmph_DetOnly = v_kmph(row);

figure; imagesc(t, v_kmph, 20*log10(abs(normalised_FFT)));
colormap('jet'); colorbar;                                                                             
xlabel('Time (s)');
ylim([2 50]);
clim([-60 0]);
ylabel('Speed (km/h)');         
axis xy; 

hold on;
p = plot(t_DetOnly, v_kmph_DetOnly, 'kx', 'markersize', 10);
p.Color = '#00841a';
x = 1; %Ask about
hold off;

%Filtering outliers from detection matrix and finding polynomial model
filtered_v_kmph_DetOnly = filloutliers(v_kmph_DetOnly, "spline");
[p,S] = polyfit(t_DetOnly, filtered_v_kmph_DetOnly, 1);

%Finding length of signal in seconds
N = length(y); % sample length
slength = N/fs;
t_samples = 0:1:round(slength);

%Evaluating polynomial at new time increments
v_poly = polyval(p, t_samples);

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
    reduced_indices = find(v_kmph <= 50 & v_kmph >= 1.6);
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