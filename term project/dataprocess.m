clear

filename = 'trace.xlsx';
sheet = 2;
xlRange = 'A1:BT3001';

subsetA = xlsread(filename,sheet,xlRange);
imagesc(subsetA(2:100,:))


plot(subsetA(2:end,6))