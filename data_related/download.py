import matplotlib.pyplot as plt
import fix_yahoo_finance as yf
import pandas as pd
d1=pd.read_csv("bse30.csv") 
comp_names=d1["Security Name"].tolist()
 
f='2017-12-18'
t='2018-09-30'
total_data=list()
for i in range(len(comp_names)):
	comp_symbol=comp_names[i].strip()+'.BO'
	print comp_symbol
	data = yf.download(comp_symbol,f,t)
	adj_close=data["Adj Close"].tolist()
	total_data.append(adj_close)
	print len(adj_close)

df=pd.DataFrame(total_data)

df = df.transpose()
df.columns = comp_names
print df
df.to_csv('final_list.csv')
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