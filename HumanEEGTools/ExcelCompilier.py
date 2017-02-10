import pandas as pd
import glob

#requires xlrd module with pandas
#requires pandas

#will merge by the specified column
##finProd = pd.merge(file, whichColumn, how='left')

#TO DO:
#add better file searching
#add usr input
#add output adjust


for f in glob.glob("ECR*"):
    print(f)

finProd = pd.DataFrame()
newDf = pd.DataFrame()
for f in glob.glob("ECR*"):
    df = pd.read_excel(open(f, 'rb'))
    finProd = finProd.append(df['CH 1, Power'], ignore_index=True)
    
finProd = finProd.transpose()
print(finProd.head())

finProd.to_csv('ECR.csv')
    #newDf = pd.merge( df['CH 1, Power'], how='right')
    #finProd = pd.conat([df['CH 1, Power'], finProd])
    #newDf = newDf.append(finProd['CH 1, Power'])
newDf = df['Frequency']
    
#=new!$B$2:$B$131,new!$C$2:$C$131,new!$D$2:$D$131,new!$E$2:$E$131,new!$F$2:$F$131,new!$G$2:$G$131,new!$H$2:$H$131,new!$I$2:$I$131, new!$J$2:$J$131,new!$K$2:$K$131, new!$L$2:$L$131,new!$M$2:$M$131,new!$N$2:$N$131
    
#finProd.describe()
#finProd.head()

#pd.to_excel('dir/newDoc.xlsx', sheet_name='Sheet1')


