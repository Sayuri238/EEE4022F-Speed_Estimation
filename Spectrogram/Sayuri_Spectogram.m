
fs = 192000; %Sampling frequency 
f1 = 20000; %Start frequency of the chirp 
f2 = 40000; %Stop frequency of the chirp
dur = 1; %Duration of chrip
ph0 = 0; %Phase not specified

%Calculate rate of hcnage of frequency
FreqChange = (f2-f1)/dur;

%Use a vector for time
t = single(0:1/fs:dur);

%Chirp signal generated
y = sin(2*pi*(f1*t + 0.5*FreqChange*t.^2));                                                              % Abdul Gaffar

%Frame Length chosen to improve frquency resolution
%Hann window chosen because Hamming window may a minor discontinuities at
%edges
FrameLen = 512;                                                                                                
ov = 0.5;    
[s, f, t] = Sayuri_Spectrogram(y, FrameLen, ov, fs);

figure; imagesc(t, f/1e3, 20*log10(abs(s)));
colormap('jet'); colorbar;                                                                             
xlabel('Number of columns');       
ylabel('Frequency (kHz)');         
axis xy; 

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
    reduced_indices = find(v_kmph <= 70 & v_kmph >= 2);
    reduced_FFT = s(reduced_indices, :);

    v_kmph_output = v_kmph(reduced_indices);

    %Finding the maximum value of FFT matrix
    max_FFT_Value = max(max((reduced_FFT)));

    %Normalising and plotting reduced matrix 
    normalised_FFT = reduced_FFT/max_FFT_Value;     
end