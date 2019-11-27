""" 
This file reads a csv file that is downloaded from the UN comtrade database, and extracts the relevant data based on 
which product is wanted based on the product code. Generates a file for each year.
Author: Steinar Laenen (steinar9@gmail.com)
"""

import csv
from tqdm import tqdm

product = "2709"
# n = 0
# with open('./type-C_r-ALL_ps-2017_freq-A_px-HS_pub-20190521_fmt-csv_ex-20190530.csv') as f:
#     n = sum(1 for line in f)
# print(n)
years = range(1988, 2018)
for year in tqdm(years):
    with open('/media/steinar/Elements/UNComtrade database/type-C_r-ALL_ps-'+ str(year) +'_freq-A_px-HS_pub.csv', mode='rb') as trade_data_file, open('./processed_data_' + str(year)+ '_' + str(product) + '.csv', mode='wb') as cleaned_file:
        reader = csv.reader(trade_data_file)

        header = reader.next()

        index_commodity_code = header.index('Commodity Code')

        writer = csv.writer(cleaned_file)
        # first write down the header
        writer.writerow(header)

        for i, row in tqdm(enumerate(reader)):
             if row[index_commodity_code] == str(product):
                 writer.writerow(row)

