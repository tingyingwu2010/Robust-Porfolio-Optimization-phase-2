# import matplotlib.pyplot as plt
# import fix_yahoo_finance as yf
import yfinance as yf
import pandas as pd
import math
import numpy as np
# d1=pd.read_csv("SnP100.csv")
# d1=pd.read_csv("SnP500.csv")
# d1=pd.read_csv("Nikkei225.csv")
# d1=pd.read_csv("JPXNikkei400.csv")
d1=pd.read_csv("ftse100.csv")
# d1=pd.read_csv("ftse350.csv")

comp_names=d1["Security Name"].tolist()
 
# f='2016-12-18'
f='2017-12-18'
t='2018-09-30'
new_array = [];
total_data=list()
for i in range(len(comp_names)):
    # comp_symbol=comp_names[i].strip()+'.BO'
    comp_symbol=comp_names[i].strip()
    # comp_symbol=comp_names[i].strip()+'.L'
    # comp_symbol=comp_names[i].strip()+'.T'
    # comp_symbol=str(comp_names[i]).strip()+'.T'
    # if (comp_symbol == 'IDEA.BO' or comp_symbol == 'CONCOR.BO' or comp_symbol == 'HDFCLIFE.BO'):
    #     continue
    print(comp_symbol)
    try:
        data = yf.download(comp_symbol,f,t)
    except:
        print("Download failed for ticker symbol: "+comp_symbol)
        continue
    # data = yf.download(comp_symbol,f,t)
    if (len(data)==0):
        print("Download failed for ticker symbol: "+comp_symbol)
        continue
    print(len(data))
    adj_close=data["Adj Close"].tolist()
    
    if True in np.isnan(adj_close):
        print("NaN values for ticker symbol: "+comp_symbol)
        continue
    # print(adj_close)
    is_price_negative=any(n < 0 for n in adj_close)
    if (is_price_negative):
        print("Negative Stock values for ticker symbol: "+comp_symbol)
        continue

   
    total_data.append(adj_close)
    print(len(adj_close))
    new_array.append(comp_symbol)
    
df=pd.DataFrame(total_data)

df = df.transpose()
df.columns = new_array
print(df)
# df.to_csv('final_snp100.csv')
# df.to_csv('final_snp500.csv')
# df.to_csv('final_nikkei225.csv')
df.to_csv('final_ftse100.csv')


#print df
#data = yf.download('INFY.BO',f,t)
#print type(data)
#print data
#adj_close=data["Adj Close"].tolist()
#a=adj_close[2]
#print(adj_close)
#print(type(a))
#print type(adj_close)
#data.Close.plot()
#plt.show()
