%% --- 4-Channel MCA Coincidence Analyzer --- %%
%%%        Author: Dr. Hamid Basiri         %%%%
%%%%%     Prof. Matsushima Laboratory     %%%%%%
%%%%%%%     The University of Tokyo     %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script analyzes Multi-Channel Analyzer (MCA) data manufactured
% by TechnoAP Co., Ltd. (APG7400A model) to identify coincidences.

% Instructions:
% 1. Record data in List mode with USB-MCA4 software.
% 2. Convert binary data to CSV via the software's file menu.
% 3. Replace FILE_NAME with your CSV file's name and path.
% 4. Adjust COIN_WINDOW to fit your timing needs for coincidence detection.
% 5. The output is displayed on the screen and recorded on a csv file
% FILE_Name_analysis.csv

% COPYRIGHT:
% © 2024 Dr. Hamid Basiri. All rights reserved. This code is for educational
% and research use only. No part may be reproduced or used for commercial
% purposes without permission from the author.

% DISCLAIMER:
% Provided "as is", without any warranty. Not liable for any damages arising
% from its use. Users assume all risks related to the performance and use of
% this software.

clc;
clear variables;
close all;

% Set the full path of the file directly
FILE_NAME = 'C:\Users\6122338134\OneDrive - The University of Tokyo\MCA_analyzer\acryl1piece_long_24hour.txt'; % Modify the path as needed
COIN_WINDOW = 1000; % in nanoseconds
OUTPUT_CSV = [strrep(FILE_NAME, '.txt', ''), '_analysis.csv']; % Output file 

% Check if the file exists
if exist(FILE_NAME, 'file') ~= 2
    error('File does not exist.');
end

fid = fopen(FILE_NAME);
if fid == -1
    error('Error opening file.');
end

% Read data from CSV
data = textscan(fid, '%f64 %d %d', 'Delimiter', ',');
fclose(fid);
timestamps = data{1};
channels = data{2};
phas = data{3};

% Initialize histograms
edges = 0:1:1024;
histCh = cell(1, 4);
for i = 1:4
    histCh{i} = zeros(length(edges), 1);
end

% For double and triple hits
histCh12 = zeros(1025, 1025);
histCh23 = zeros(1025, 1025);
histCh31 = zeros(1025, 1025);
histCh14 = zeros(1025, 1025); 
histCh34 = zeros(1025, 1025); 

histCh123 = zeros(1025, 1025, 1025);
histCh134 = zeros(1025, 1025, 1025);  
histCh124 = zeros(1025, 1025, 1025);  
histCh234 = zeros(1025, 1025, 1025);  

quadraHitsCount = 0;
noOfNHits = zeros(1, 5); 


% Calculate total measurement time in minutes and hours
totalTimeInMinutes = (timestamps(end) - timestamps(1)) * 1e-9 / 60;
totalTimeInHours = totalTimeInMinutes / 60;

% Determine number of complete hours in dataset
numHours = ceil(totalTimeInHours);  % Calculate the number of hours covered

% Initialize arrays to count hits per hour for each type
hourlyCountsSingle = zeros(1, numHours);
hourlyCountsDouble = zeros(1, numHours);
hourlyCountsTriple = zeros(1, numHours);
hourlyCountsQuadra = zeros(1, numHours);

% Process events
prevTime = -COIN_WINDOW;
nHits = 0;
nEvent = 0;
pha = zeros(1, 4);
time = zeros(1, 4);

for i = 1:length(timestamps)
    currTime = timestamps(i);
    nCh = channels(i); % Use direct channel number from CSV
    tempPha = min(phas(i), 1024);
    hourIndex = floor(((currTime - timestamps(1)) * 1e-9) / 3600) + 1; % Hour index

    if currTime - prevTime < COIN_WINDOW && prevTime ~= -COIN_WINDOW
        nHits = nHits + 1;
        pha(nCh) = tempPha;
        time(nCh) = currTime;
        histCh{nCh}(tempPha + 1) = histCh{nCh}(tempPha + 1) + 1;

        % Check and record new histograms
        if nHits >= 2
            if pha(1) > 0 && pha(4) > 0
                histCh14(pha(1) + 1, pha(4) + 1) = histCh14(pha(1) + 1, pha(4) + 1) + 1;
            end
            if pha(3) > 0 && pha(4) > 0
                histCh34(pha(3) + 1, pha(4) + 1) = histCh34(pha(3) + 1, pha(4) + 1) + 1;
            end
        end
    else
        % Process triple coincidences before resetting
        if nHits == 3
            if pha(1) > 0 && pha(3) > 0 && pha(4) > 0
                histCh134(pha(1) + 1, pha(3) + 1, pha(4) + 1) = histCh134(pha(1) + 1, pha(3) + 1, pha(4) + 1) + 1;
            end
            if pha(1) > 0 && pha(2) > 0 && pha(4) > 0
                histCh124(pha(1) + 1, pha(2) + 1, pha(4) + 1) = histCh124(pha(1) + 1, pha(2) + 1, pha(4) + 1) + 1;
            end
            if pha(2) > 0 && pha(3) > 0 && pha(4) > 0
                histCh234(pha(2) + 1, pha(3) + 1, pha(4) + 1) = histCh234(pha(2) + 1, pha(3) + 1, pha(4) + 1) + 1;
            end
        end

        if nHits == 4
            quadraHitsCount = quadraHitsCount + 1; % Correctly increment quadra hit count
            hourlyCountsQuadra(hourIndex) = hourlyCountsQuadra(hourIndex) + 1;
        end
        % Process the previous event if there was one
        if nHits > 0
            nHits = min(nHits, length(noOfNHits));  % Cap nHits at the size of the array
            noOfNHits(nHits) = noOfNHits(nHits) + 1;
            nEvent = nEvent + 1;
            if hourIndex <= numHours
                switch nHits
                    case 1
                        hourlyCountsSingle(hourIndex) = hourlyCountsSingle(hourIndex) + 1;
                    case 2
                        hourlyCountsDouble(hourIndex) = hourlyCountsDouble(hourIndex) + 1;
                    case 3
                        hourlyCountsTriple(hourIndex) = hourlyCountsTriple(hourIndex) + 1;
                end
            end
        end
        % Reset for the new event
        nHits = 1;
        pha = zeros(1, 4);
        time = zeros(1, 4);
        pha(nCh) = tempPha;
        time(nCh) = currTime;
        prevTime = currTime;
        histCh{nCh}(tempPha + 1) = histCh{nCh}(tempPha + 1) + 1;
    end
end


% Finalize the last event
if nHits > 0
    nHits = min(nHits, length(noOfNHits));  % Cap nHits at the size of the array
    noOfNHits(nHits) = noOfNHits(nHits) + 1;
end

% Calculate cpm and cph for each event type
cpmSingle = noOfNHits(1) / totalTimeInMinutes;
cpmDouble = noOfNHits(2) / totalTimeInMinutes;
cpmTriple = noOfNHits(3) / totalTimeInMinutes;
cpmQuadra = quadraHitsCount / totalTimeInMinutes;

cphSingle = noOfNHits(1) / totalTimeInHours;
cphDouble = noOfNHits(2) / totalTimeInHours;
cphTriple = noOfNHits(3) / totalTimeInHours;
cphQuadra = quadraHitsCount / totalTimeInHours;

% Calculate the standard deviation for each type of hit event per hour
stdSingle = std(hourlyCountsSingle);
stdDouble = std(hourlyCountsDouble);
stdTriple = std(hourlyCountsTriple);
stdQuadra = std(hourlyCountsQuadra);

% Uncertainty calculation using sqrt(count)
uncertaintySingle = sqrt(cphSingle);
uncertaintyDouble = sqrt(cphDouble);
uncertaintyTriple = sqrt(cphTriple);
uncertaintyQuadra = sqrt(cphQuadra);

% Prepare data for output table including uncertainties
summaryTable = table({'Single Hit'; 'Double Hit'; 'Triple Hit'; 'Quadruple Hit'}, ...
                     [noOfNHits(1); noOfNHits(2); noOfNHits(3); quadraHitsCount], ...
                     [cpmSingle; cpmDouble; cpmTriple; cpmQuadra], ...
                     [cphSingle; cphDouble; cphTriple; cphQuadra], ...
                     [stdSingle; stdDouble; stdTriple; stdQuadra], ...
                     [uncertaintySingle; uncertaintyDouble; uncertaintyTriple; uncertaintyQuadra], ...
                     'VariableNames', {'Event_Type', 'Total_Hits', 'Counts_Per_Minute', 'Counts_Per_Hour', 'Standard_Deviation_Per_Hour', 'Uncertainty_Per_Hour'});

writetable(summaryTable, OUTPUT_CSV);


disp(['Number of events: ', num2str(nEvent)]);
disp(summaryTable);

%% Drawing histograms for each channel
figure(WindowState="maximized");
for k = 1:4
    subplot(2, 2, k);
    bar(histCh{k});
    maxVal = max(histCh{k});
    roundedMax = max(100 * ceil(maxVal / 100), 100);  % Ensure minimum range of 100
    ylim([0 roundedMax])
    title(['Ch', num2str(k), ' Pulse Height']);
end

%% Determine the maximum value from all channels for unified plotting
allMaxVals = cellfun(@max, histCh);
maxOverall = max(allMaxVals);
roundedOverallMax = max(100 * ceil(maxOverall / 100), 100);

% Drawing all channels in one plot with lines
figure(WindowState="maximized");
subplot(2,1,1); % Subplot for raw data
hold on;
plot(histCh{1}, 'k', 'LineWidth', 2);  
plot(histCh{2}, 'r', 'LineWidth', 2);  
plot(histCh{3}, 'g', 'LineWidth', 2);  
plot(histCh{4}, 'b', 'LineWidth', 2);  
hold off;
legend('Channel 1', 'Channel 2', 'Channel 3', 'Channel 4');
title('Raw Pulse Heights for All Channels');
xlabel('Pulse Height');
ylabel('Count');
xlim([0 1024]);
ylim([0 roundedOverallMax]);  

% Subplot for smoothed data
subplot(2,1,2);
hold on;
smoothingSpan = 10;  % Adjust the smoothing window size as needed
plot(smooth(histCh{1}, smoothingSpan), 'k', 'LineWidth', 2);  
plot(smooth(histCh{2}, smoothingSpan), 'r', 'LineWidth', 2);  
plot(smooth(histCh{3}, smoothingSpan), 'g', 'LineWidth', 2);  
plot(smooth(histCh{4}, smoothingSpan), 'b', 'LineWidth', 2);  
hold off;
legend('Channel 1', 'Channel 2', 'Channel 3', 'Channel 4');
title('Smoothed Pulse Heights for All Channels');
xlabel('Pulse Height');
ylabel('Count');
xlim([0 1024]);
ylim([0 roundedOverallMax]);  
