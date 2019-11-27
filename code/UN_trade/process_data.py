import csv

with open('./trade_data.csv', mode='rb') as trade_data_file, open('./cleaned_data.csv', mode='wb') as cleaned_file:
    writer = csv.writer(cleaned_file)
    first = True
    
    for row in csv.reader(trade_data_file):
        if 'No data' not in row[0] and row[0] != 'Classification' or first == True:
            writer.writerow(row)
            first = False


