close all;
clear; clc;
%% Setting up Directory
% Define paths for the project, data, and results directories
parentPath = 'C:\Users\brth229\OneDrive - University of Kentucky\Research Projects\ECGDetection\'; %Path of the project folder
dataPath = [parentPath, 'data\']; % (https://www.kaggle.com/dsv/6114424)
resultsPath = [parentPath, 'results\'];
%% Creating a sub direcotry in for MLII results 
resultsPathMLII = [resultsPath,'MLII\']; % Rename if needed for different lead types.
[~,~,~] = mkdir(resultsPathMLII); % Ensure the results directory exists
%% Extracting all data files from dataPath
% Get all .csv files in the directory and its subdirectories
csvFiles = dir(fullfile(dataPath, '**/*.csv'));
% Sort files based on numbers in filenames
[~, sortIdx] = sort(cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), {csvFiles.name}));
% Extract the file names from the list of file information structures
% Apply sorting
csvFiles = csvFiles(sortIdx);
recordingFileNames = {csvFiles.name};
%% Global Variables
% Sampling frequency of ECG signals (from the dataset documentation)
fs = 360; 
% Bandpass filter cutoff frequncies
lowCutoff = 8; highCutoff = 20;
%% Creating a cell array for results
% Define the column names for the results table
ColumnNames = {'Data #', 'Total Beats','Detected R Peaks', 'Detected QRS Beats','QRS Detection Rate',...
    'Number of WCT', 'Number of NCT', 'Average QRS Width', 'Standard Deviation'};
% Initialize a cell array to store the results
resultCell = cell(1,length(ColumnNames));
resultCell(1,:) = ColumnNames;
%% Total beats from the dataset excluding 102 and 104 (46 files)
% Derived from:
% https://www.physionet.org/files/mitdb/1.0.0/mitdbdir/records.htm
totalBeats = [2273,1865,2084,2572,2027,2137,1774,2532,2124,2539,1795,1879,...
    1953,2412,1535,2288,1987,1863,2476,1518,1619,2601,2000,2136,2980,2656,2332,2955,...
    3005,2650,2748,3251,2262,3363,2208,2287,2048,2427,2483,2605,2053,2256,1573,1780,3079,2753];
%% Loop over all recordings
%Reads patient data from our 3D Array and plots respective graphs
for irecording = 1:length(recordingFileNames)
    %% Extracting ECG data
    % Get the current file name and corresponding number
    currentfileName = recordingFileNames{irecording};
    fileNumber = split(currentfileName,'.');
    fileNumber = fileNumber{1};
    % Load the ECG data, skipping the first row as it stores the headers
    % for column names.
    completeEKGData = readmatrix(fullfile(dataPath,currentfileName), 'Range', 2);
    disp(['<strong>Data ID: ',currentfileName, ' is extracted.</strong>'])
    % Extract variable names (headers) to identify the columns
    headers = string(readtable(fullfile(dataPath,currentfileName)).Properties.VariableNames);
    timeIndex = find(headers == 'time_ms'); % Column index for time
    MLIIIndex = find(headers=='MLII'); % Column index for MLII lead
    V5Index = find(headers == 'V5'); % Column index for V5 lead (if available)
    if isempty(V5Index)
        V5Index=4; % Default to 4 if V5 lead is not present
    end
    %% Extracting relevant leads and time data
    time = completeEKGData(:,timeIndex)'; % Time data in milliseconds
    L = length(time); % Length of the time series
    MLII = completeEKGData(:,MLIIIndex)'; % MLII lead data
    V5 = completeEKGData(:,V5Index)'; % V5 lead data (if available)
    %% Filtering the ECG signal
    [filteredMLII,filterMLII] = bandpass(MLII, [lowCutoff highCutoff], fs,'ImpulseResponse','fir');
    [filteredV5, filterV5] = bandpass(V5, [lowCutoff highCutoff], fs,'ImpulseResponse','fir');
    %% Unfiltered MLII Visualization
    % Configure plot properties for time and frequency-domain visualizations
    commonLineWidth = 1;
    commonFontSize = 8;
    figureWidth = 3.5;
    eegHeight = 1.5;
    commonXLimTime = [90,100];    commonXLimFFT = [0,180];
    commonYLimTime = [-1,1];    commonYLimFFT = [0,0.01];
    % Compute the FFT of the unfiltered MLII signal
    fftMLII = fft(MLII);
    P1 = abs(fftMLII/L);
    P2 = P1(1:floor(L/2)+1);
    P2(2:end-1) = 2*P2(2:end-1);
    f = fs*(0:(L/2))/L;

    % Save time-domain and FFT visualizations for the unfiltered MLII lead
    MLIITimeFigure = figure('Visible','off');
    set(MLIITimeFigure, 'PaperUnits', 'inches', 'PaperPosition', [0 0 figureWidth eegHeight]);
    plot(time./1000, MLII,'LineWidth',commonLineWidth)
%     title(['Recording ' num2str(currentfileName) ' MLII Lead Data'])
    xlabel("Time [s]",'FontWeight','bold')
    ylabel("Amplitude [\muV]",'FontWeight','bold');yticks([-2:2]);yticklabels([-2:2])
    set(gca,'YLim',commonYLimTime,'XLim',commonXLimTime,'FontSize',commonFontSize)
    fileName = ['Fig_' fileNumber '_MLII_Time' ];
    saveas(gcf, fullfile(resultsPathMLII, [fileName '.fig']));
    print(gcf, fullfile(resultsPathMLII, [fileName '.tif']), '-dtiff', '-r300');


    MLIIFreqFigure = figure('Visible','off');
    set(MLIIFreqFigure, 'PaperUnits', 'inches', 'PaperPosition', [0 0 figureWidth eegHeight]);
    plot(f, P2,'LineWidth',1)
%     title(['Patient ' num2str(currentfileName) ' MLII FFT'])
    xlabel("Frequency [Hz]",'FontWeight','bold')
    ylabel("Magnitude",'FontWeight','bold')
    xlim([0, 180]); xticks([0:30:180]);yticklabels([0:30:180])
    ylim([0, 0.01])
    yticks([0,0.01]);yticklabels([0,0.01])
    set(gca,'YLim',commonYLimFFT,'XLim',commonXLimFFT,'FontSize',commonFontSize)
    fileName = ['Fig_' fileNumber '_MLII_FFT' ];
    saveas(gcf, fullfile(resultsPathMLII, [fileName '.fig']));
    print(gcf, fullfile(resultsPathMLII, [fileName '.tif']), '-dtiff', '-r300');

    close all;
    %% Filtered MLII Visualization
    % Compute the FFT of the filtered MLII signal
    fftBPMLII = fft(filteredMLII);
    BP1 = abs(fftBPMLII/L);
    BP2 = BP1(1:floor(L/2)+1);
    BP2(2:end-1) = 2*BP2(2:end-1);

    % Save time-domain and FFT visualizations for the filtered MLII lead
    MLIITimeFigure = figure('Visible','off');
    set(MLIITimeFigure, 'PaperUnits', 'inches', 'PaperPosition', [0 0 figureWidth eegHeight]);
    plot(time./1000, filteredMLII,'LineWidth',commonLineWidth)
%     title(['Recording ' num2str(currentfileName) ' MLII Lead Data'])
    xlabel("Time [s]",'FontWeight','bold')
    ylabel("Amplitude [\muV]",'FontWeight','bold');yticks([-2:2]);yticklabels([-2:2])
    set(gca,'YLim',commonYLimTime,'XLim',commonXLimTime,'FontSize',commonFontSize)
    fileName = ['Fig_' fileNumber '_FilteredMLII_Time'];
    saveas(gcf, fullfile(resultsPathMLII, [fileName '.fig']));
    print(gcf, fullfile(resultsPathMLII, [fileName '.tif']), '-dtiff', '-r300');


    MLIIFreqFigure = figure('Visible','off');
    set(MLIIFreqFigure, 'PaperUnits', 'inches', 'PaperPosition', [0 0 figureWidth eegHeight]);
    plot(f, BP2,'LineWidth',commonLineWidth)
    title(['Patient ' num2str(currentfileName) ' MLII FFT'])
    xlabel("Frequency [Hz]",'FontWeight','bold')
    ylabel("Magnitude [\muV]",'FontWeight','bold')
    xlim([0, 180])
    ylim([0, 0.01])
    set(gca,'YLim',commonYLimFFT,'XLim',commonXLimFFT,'FontSize',commonFontSize)
    fileName = ['Fig_' fileNumber '_FilteredMLII_FFT' ];
    saveas(gcf, fullfile(resultsPathMLII, [fileName '.fig']));
    print(gcf, fullfile(resultsPathMLII, [fileName '.tif']), '-dtiff', '-r300');
    close all;

    % V5 channel
    V5TimeFigure = figure('Visible','off');
    set(V5TimeFigure, 'PaperUnits', 'inches', 'PaperPosition', [0 0 figureWidth eegHeight]);
    plot(time./1000, filteredV5,'LineWidth',commonLineWidth)
    %     title(['Recording ' num2str(currentfileName) ' MLII Lead Data'])
    xlabel("Time [s]",'FontWeight','bold')
    ylabel("Amplitude [\muV]",'FontWeight','bold');yticks([-2:2]);yticklabels([-2:2])
    set(gca,'YLim',commonYLimTime,'XLim',commonXLimTime,'FontSize',commonFontSize)
    fileName = ['Fig_' fileNumber '_FilteredV5_Time'];
    saveas(gcf, fullfile(resultsPathMLII, [fileName '.fig']));
    print(gcf, fullfile(resultsPathMLII, [fileName '.tif']), '-dtiff', '-r300');
    %% Filtered MLII Visualization for 108 and 116
    if str2double(fileNumber) == 108
        % Save time-domain and FFT visualizations for the filtered MLII lead
        MLIITimeFigure = figure('Visible','on');
        set(MLIITimeFigure, 'PaperUnits', 'inches', 'PaperPosition', [0 0 figureWidth eegHeight]);
        plot(time./1000, filteredMLII,'LineWidth',commonLineWidth)
        %     title(['Recording ' num2str(currentfileName) ' MLII Lead Data'])
        xlabel("Time [s]",'FontWeight','bold')
        ylabel("Amplitude [\muV]",'FontWeight','bold');yticks([-2:2]);yticklabels([-2:2])
        set(gca,'YLim',commonYLimTime,'XLim',[445,475],'FontSize',commonFontSize)
        yline(0.2,'LineWidth',commonLineWidth)
        fileName = ['Fig_' fileNumber '_FilteredMLII_Time'];
        saveas(gcf, fullfile(resultsPathMLII, [fileName '.fig']));
        print(gcf, fullfile(resultsPathMLII, [fileName '.tif']), '-dtiff', '-r300');
    elseif str2double(fileNumber) == 116
        % Save time-domain and FFT visualizations for the filtered MLII lead
        MLIITimeFigure = figure('Visible','on');
        set(MLIITimeFigure, 'PaperUnits', 'inches', 'PaperPosition', [0 0 figureWidth eegHeight]);
        plot(time./1000, filteredMLII,'LineWidth',commonLineWidth)
        %     title(['Recording ' num2str(currentfileName) ' MLII Lead Data'])
        xlabel("Time [s]",'FontWeight','bold')
        ylabel("Amplitude [\muV]",'FontWeight','bold');yticks([-2:2]);yticklabels([-2:2])
        set(gca,'YLim',commonYLimTime,'XLim',[440,450],'FontSize',commonFontSize)
        yline(0.2,'LineWidth',commonLineWidth)
        fileName = ['Fig_' fileNumber '_FilteredMLII_Time'];
        saveas(gcf, fullfile(resultsPathMLII, [fileName '.fig']));
        print(gcf, fullfile(resultsPathMLII, [fileName '.tif']), '-dtiff', '-r300');

    end


    %% QRS Detection
    % Detect R-peaks and compute QRS complex features (QRS width, Q, and S points)
    % Set parameters for peak detection
    peakThreshold = 0.2;  % Minimum height of a peak for it to be considered an R-peak
    peakDistance = round(0.2*fs); % Minimum distance between consecutive R-peaks in samples (200 ms)

    % Detect R-peaks in the filtered MLII signal
    [rPeaks, rLocs] = findpeaks(filteredMLII, 'MinPeakHeight', peakThreshold, 'MinPeakDistance',peakDistance);


    % Initialize arrays to store Q, S points and their amplitudes
    qIndices = zeros(length(rPeaks), 1);
    sIndices = zeros(length(rPeaks), 1);
    qAmps = zeros(length(rPeaks), 1);
    sAmps = zeros(length(rPeaks), 1);

    % Define QRS window duration based on a reference from the literature
    % Reference: QRS interval for a healthy adult is 100 ± 20 ms
    % Source: Mohamed Elgendi, Mirjam Jonkman, and Friso DeBoer
    qrsWindow=round(0.120*fs); %120 ms Window
    qrsWidths = zeros(length(rPeaks), 1);  % Array to store QRS widths

    % Initialize variables to find local maxima around Q and S points
    maxBeforeQ = zeros(length(rPeaks), 1);
    maxAfterS = zeros(length(rPeaks), 1);
    maxBeforeQIndices = zeros(length(rPeaks), 1);
    maxAfterSIndices = zeros(length(rPeaks), 1);
    qrsWidthsMaxima = zeros(length(rPeaks), 1);

    totalQRSCounter=0; % Counter for total detected QRS complexes
    % Loop through each detected R-peak to identify Q, S points and calculate QRS widths
    for j = 1:length(rPeaks)
        % Define the index of the current R-peak
        rPeakIndex = rLocs(j); 

        % Define the search window for Q point (before R-peak)
        startQ = max(1, rPeakIndex - qrsWindow);  % Ensure within bounds
        endQ = rPeakIndex - 1;  % End just before R-peak

        % Set the window for S (after R-peak)
        startS = rPeakIndex + 1;  % Start just after R-peak
        endS = min(length(filteredMLII), rPeakIndex + qrsWindow);  % Ensure within bounds

        % Find Q point (minimum amplitude in the Q window)
        [qAmps(j), qIndices(j)] = min(filteredMLII(startQ:endQ));
        % Find S point (minimum amplitude in the S window)
        [sAmps(j), sIndices(j)] = min(filteredMLII(startS:endS));

        % Map the Q and S indices back to the original signal
        qIndices(j) = qIndices(j) + startQ - 1; % Adjust for the startQ offset
        sIndices(j) = sIndices(j) + startS - 1; % Adjust for the startS offset
        % Calculate QRS width as the time difference between S and Q points
        qrsWidths(j) = (time(sIndices(j)) - time(qIndices(j)));  % QRS width in milliseconds


        % Find the local maxima after the S minima
        for k = qIndices(j):-1:2   % Iterate backwards from the Q point
            if filteredMLII(k - 1) < filteredMLII(k) && filteredMLII(k) > filteredMLII(k + 1)
                maxBeforeQIndices(j) = k; % Index of local maximum
                maxBeforeQ(j) = filteredMLII(k); % Amplitude of local maximum
                break;
            end
        end

        % Find the local maxima after the S minima
        for k = sIndices(j):length(filteredMLII) - 1  % Iterate forwards from the S point
            if filteredMLII(k - 1) < filteredMLII(k) && filteredMLII(k) > filteredMLII(k + 1)
                maxAfterSIndices(j) = k; % Index of local maximum
                maxAfterS(j) = filteredMLII(k); % Amplitude of local maximum
                break;
            end
        end



        % Calculate QRS width using the local maxima around the QRS complex
        if maxBeforeQIndices(j) > 0 && maxAfterSIndices(j) > 0
            qrsWidthsMaxima(j) = time(maxAfterSIndices(j)) - time(maxBeforeQIndices(j));  % Width in milliseconds
            totalQRSCounter = totalQRSCounter +1; % Increment the QRS counter
        else
            qrsWidthsMaxima(j) = NaN;  %% Assign NaN if maxima are not found
        end


    end

    %% Compute Average and STD and if Exceeds 120ms (WCT)

    % Calculate the mean of QRS widths
    % 'omitnan' ensures that NaN values (from undetected QRS complexes) are ignored
    meanQRSWidthMaxima = mean(qrsWidthsMaxima, 'omitnan');  % Mean, ignoring NaN values
    stdQRSWidthMaxima = std(qrsWidthsMaxima, 'omitnan');    % Standard deviation, ignoring NaN values

    % Find the indices of QRS complexes whose widths exceed 120 ms
    % A QRS width >120 ms is indicative of Wide Complex Tachycardia (WCT)
    exceedIndices = find(qrsWidthsMaxima > 120);
    % Initialize a counter to track the number of Wide Complex Tachycardia cases
    wctCounter=0;
    % Check if there are any QRS complexes exceeding 120 ms
    if isempty(exceedIndices)
        % Display a message if no QRS complexes exceed the threshold
        disp('No QRS widths exceed 120 ms.');
    else
        % disp('QRS widths exceeding 120 ms are found at the following indices:');
        for j = 1:length(exceedIndices)
            wctCounter = wctCounter + 1;
        end
    end
    %% Compute Detection Accuracy
    detectionAccuracy = round((1 - abs(totalBeats(1,irecording)-totalQRSCounter)/totalBeats(1,irecording))*100,2);

    %% Display and Storage of Summary
    % Display results
    disp(['Total Beats: ' num2str(totalQRSCounter)]);
    disp(['Average QRS Width: ', num2str(meanQRSWidthMaxima), ' ms']);
    disp(['Standard Deviation: ', num2str(stdQRSWidthMaxima), ' ms']);
    disp(['Detection Accuracy: ', num2str(detectionAccuracy), '%']);
    disp(['Number of WCT: ', num2str(wctCounter), ' out of ', num2str(totalQRSCounter), ' total QRS complexes']);
    disp(['Number of NCT: ', num2str(totalQRSCounter-wctCounter), ' out of ', num2str(totalQRSCounter), ' total QRS complexes']);

    % Save results in resultCell in every loop
    resultCell(1+irecording,:) = {fileNumber, totalBeats(1,irecording) ,length(rPeaks), totalQRSCounter,detectionAccuracy ,wctCounter, totalQRSCounter-wctCounter, meanQRSWidthMaxima, stdQRSWidthMaxima} ;

    %% Figure for detected QRS
    detectedECGFigure = figure('Visible','off');
    set(detectedECGFigure, 'PaperUnits', 'inches', 'PaperPosition', [0 0 3.5 1.5]);
    plot(time./1000, filteredMLII,'LineWidth',1);
    hold on;
    plot(time(rLocs)./1000, rPeaks, 'ro','MarkerSize',4);  % R peaks in red
    plot(time(qIndices)./1000, filteredMLII(qIndices), 'go','MarkerSize',4);  % Q points in green
    plot(time(sIndices)./1000, filteredMLII(sIndices), 'bo','MarkerSize',4);  % S points in blue
    xLimDetected = [90 95];
    set(gca,'YLim',commonYLimTime,'XLim',xLimDetected)

    % Add brackets and annotate the QRS widths
    for j = 1:length(qIndices)
        if ~isnan(qrsWidthsMaxima(j)) % Skip if QRS width is not a number
            % Set the x-coordinates for the bracket
            maxBeforeQTime = time(maxBeforeQIndices(j));
            maxAfterSTime = time(maxAfterSIndices(j));
            midTime = (maxBeforeQTime + maxAfterSTime) / 2;  % Midpoint for text

            % Set the y-coordinate for the bracket height
            bracketHeight = max(rPeaks) * 1; % 1 for 100; 0.8 for 101

            if (maxBeforeQTime>=xLimDetected(1)*1000 && maxAfterSTime<=xLimDetected(2)*1000) % lines restricted to only xLims
                % Draw the bracket around the QRS complex
                line([maxBeforeQTime maxBeforeQTime]./1000, [maxBeforeQ(j) bracketHeight], 'Color', 'black','LineWidth',1);  % Left vertical line
                line([maxAfterSTime maxAfterSTime]./1000, [maxAfterS(j) bracketHeight], 'Color', 'black','LineWidth',1);  % Right vertical line
                line([maxBeforeQTime maxAfterSTime]./1000, [bracketHeight bracketHeight], 'Color', 'black','LineWidth',1);  % Top horizontal line

                % Width displayed
%                 text(midTime./1000, bracketHeight * 1.15, sprintf('%.1f ms', qrsWidthsMaxima(j)),'HorizontalAlignment', 'center', 'Color', 'black','FontSize',5);
            end
        end
    end

    xlabel("Time [s]",'FontWeight','bold','FontSize',9)
    ylabel("Amplitude [\muV]",'FontWeight','bold','FontSize',9)
    yticks([-2:2]);yticklabels([-2:2])
    lgd = legend({'ECG Signal', 'R', 'Q', 'S'},'Location','southeast','FontSize',6, 'Orientation','horizontal');
    % Adjust the gap between the marker and text
    set(lgd, 'ItemTokenSize', [5, 1]); % Smaller width (first value) reduces the gap
    set(gca,'FontSize',commonFontSize);
    fileName = ['Fig_' fileNumber '_QRSDetected_Time' ];
    saveas(gcf, fullfile(resultsPathMLII, [fileName '.fig']));
    print(gcf, fullfile(resultsPathMLII, [fileName '.tif']), '-dtiff', '-r300');
    close all;

end
%% Write the table to an Excel file
writecell(resultCell, fullfile(resultsPath,'DataResults.xlsx'));
disp(['The results are saved in:',newline, resultsPath])
disp('The execution of the script has completed.')
