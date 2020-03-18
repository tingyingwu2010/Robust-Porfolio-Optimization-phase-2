import matplotlib.pyplot as plt
import fix_yahoo_finance as yf
import pandas as pd
d1=pd.read_csv("BSE100.csv") 
comp_names=d1["Security Name"].tolist()
 
f='2016-12-18'
t='2018-09-30'
new_array = [];
total_data=list()
for i in range(len(comp_names)):
	comp_symbol=comp_names[i].strip()+'.BO'
	if (comp_symbol == 'IDEA.BO' or comp_symbol == 'CONCOR.BO' or comp_symbol == 'HDFCLIFE.BO'):
		continue
	print comp_symbol
	data = yf.download(comp_symbol,f,t)
	adj_close=data["Adj Close"].tolist()
	total_data.append(adj_close)
	print len(adj_close)
	new_array.append(comp_symbol)

df=pd.DataFrame(total_data)

df = df.transpose()
df.columns = new_array
print df
df.to_csv('final_list100.csv')
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