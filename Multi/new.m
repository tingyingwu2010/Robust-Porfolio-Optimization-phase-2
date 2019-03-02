% Read the csv file (Change it to "final_list100.csv" for BSE100)
toread = readtable('./data_related/final_list100.csv');
% table = readtable('./data_related/final_list100.csv');
stock_prices = toread{:,2:end};

fmin

