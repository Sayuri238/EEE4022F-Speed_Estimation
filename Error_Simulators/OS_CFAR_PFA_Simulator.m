clear all;
close all;

NumRows = 2000;
NumCols = 2000;
y_complex_matrix = randn(NumRows, NumCols) + 1i*randn(NumRows, NumCols); %Complex data (I and Q values)

PFA = 10^(-3); 
RefWindow = 32;
guardCells = 2; %Specifies total number of guard cells
RefWindowHalf = RefWindow/2; %Finding number of reference cells on ecah side of CUT
GuardCellsHalf = guardCells/2; %Number of guard cells on each side of CUT
CUT_index = RefWindow/2 + GuardCellsHalf + 1; %Finding index of CUT
sliderLen = RefWindow + guardCells + 1; %Creating a new 'window' to loop signal

% Adding the same number of cells as the slider length at the end of the signal 
% matrix to wrap around signal when looping
y_complexExtended = [y_complex_matrix; y_complex_matrix([1:sliderLen], :)];

%Power law detector
DataAfterPowerLawDetector_matrix = abs(y_complexExtended).^2;

%Finding length of extended matrix
extended_NumRows = NumRows + sliderLen;

detectMatrix = [];
totalNoDetections = 0;

%Looping through columns of matrix to apply OS-CFAR
for r = 1:NumCols
   DataAfterPowerLawDetector_column = DataAfterPowerLawDetector_matrix(:, r);
   [detectCol, numberOfDetections] = Sayuri_OSCFAR(RefWindow, guardCells, sliderLen, PFA, DataAfterPowerLawDetector_column, extended_NumRows, CUT_index);
   detectMatrix(:, r) = transpose(detectCol);
   totalNoDetections = totalNoDetections + numberOfDetections;
end

% Replacing the first rows that had not been used as CUT with the last rows
% of detectMatrix (so that detect points actually match up) 
detectMatrix([1:(RefWindowHalf+1)], :) = detectMatrix([(NumRows+1):(NumRows+RefWindowHalf+1)], :); 

% Replacing the last rows of extended matrix with zeros because they have
% been copied to beginning of matrix. Otherwise indices higher than NumRows
% are detected which are longer than v_kmph (access error caused)
detectMatrix([NumRows+1:(NumRows+RefWindowHalf+1)], :) = zeros((RefWindowHalf+1), NumCols);

%Finding simulation's PFA 
PFA_Simulation = totalNoDetections/(NumRows*NumCols) 

%Calculating PFA error between simualtion and defined PFA
PFA_error = ((PFA - PFA_Simulation)/PFA) * 100 % ((0.001 - 0.001)/0.001) * 100
[row, col] = find(detectMatrix == 1);

%Plotting matrix of detections
figure; imagesc(detectMatrix);
colormap('jet'); colorbar; 
hold on; 
p = plot(row, col, 'kx', 'markersize', 10);
p.Color = '#00841a';
axis xy;

function [detectCol, numberOfDetections] = Sayuri_OSCFAR(RefWindow, guardCells, sliderLen, PFA, DataAfterPowerLawDetector_column, extended_NumRows, CUT_index)
    RefWindowHalf = RefWindow/2; %Finding number of reference cells on each side of CUT
    i = 1; %index to find start of next frame
    indexBeforeExcludedCells = RefWindowHalf; %To find index before
    %guard cells and CUT
    indexAfterExcludedCells = RefWindowHalf + guardCells + 1 + 1; % To find index
    %after guard cells and CUT
    sliderEnd = sliderLen; %Last index of slider

    alphaCA = 6.0860; %Scaling factor
    gCA = [];
    TCA = []; %CA-CFAR threshold
    refList = [];
    sorted_refList = [];
    NumRows = extended_NumRows - sliderLen;
    k = floor((3*RefWindow)/4);

    while sliderEnd < extended_NumRows 
        % Defines a list without CUT or guard cells and adds it as a column
        refList(:, i) = [DataAfterPowerLawDetector_column(i:indexBeforeExcludedCells);  ...
        DataAfterPowerLawDetector_column(indexAfterExcludedCells:sliderEnd)];
        % Sorts the new column in ascending order
        sorted_refList(:, i) = sort(refList(:, i));

        %The kth value is chosen
        gCA(:, i) = sorted_refList(k, i);
       
        % Threshold is set by multiplying the scaling factor and noise
        % estimate
        TCA(:, i) = gCA(:, i)*alphaCA; 
        
        i = i + 1; %Update start of next frame length    
        indexBeforeExcludedCells = indexBeforeExcludedCells + 1;       
        indexAfterExcludedCells = indexAfterExcludedCells + 1;      
        sliderEnd = sliderEnd + 1; %Update end of next frame length     
    end
   
    numberOfDetections = 0;
    counter = 1;
    detect_indices = [];
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