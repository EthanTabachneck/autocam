import numpy as np
import pandas as pd

while(True):
    n=input('Enter amount of coefficients\n')
    if n == 'q':
        quit()
    n=int(n)
    if n<4:
        print("The amount of coefficients you have asked for is too low")
    elif n>8:
        print("The amount of coefficients you have asked for is too high")
    else:
        break

while(True):
    filename=input('Enter Path for the restrictions on the motion\n')
    if filename == 'q':
        quit()
    #Try to read the file
    try:
        res=pd.read_csv(filename)
        break
    except:
        print("The file you specified could not be found, try again")

#Convert datafram to float
for i in range(len(res)):
    res.loc[i,'w']=float(res.loc[i,'w'])
    res.loc[i,'theta'] = float(res.loc[i,'theta'])
    res.loc[i,'s'] = float(res.loc[i,'s'])
    res.loc[i,'v'] = float(res.loc[i,'v'])
    res.loc[i,'a'] = float(res.loc[i,'a'])
    res.loc[i,'j'] = float(res.loc[i,'j'])
#Convert degrees to radians for no other reason than that they are better
res['theta']=np.radians(res['theta'])
res['w']=np.radians(res['w'])

w=res['w'][0]

#creat empty coefficient array
coeffs=[]

#Matrix used for the linear algebra
M=np.array([[24,60,120,210],[12,20,30,42],[4,5,6,7],[1,1,1,1]])
#resize for the problem at hand
Mres=[]
Mfli=[]
#Make it the right size
for i in range(n-4):
    Mfli.append([])
    Mres.append([])
    for j in range(n-4):
        Mres[i].append([])
        Mfli[i].append([])
for i in range(n-4):
    for j in range(n-4-i):
        Mfli[i][j+i]=M[3-i][j+i]
        Mfli[j+i][i]=M[3-j-i][i]
        Mres[i]=Mfli[len(Mfli)-i-1]
#convert to numpy array so it can be inverted
Mres=np.array(Mres)
#Invert it
try:
    Minv=np.linalg.inv(Mres)
except:
    #This is the case where 5 coefficients are being found. Can't inverse a
    #single value
    Minv=1/Mres[0][0]

#Making a set of arrays to contain all the coefficient information
for i in range(len(res)):
    coeffs.append([])
    for j in range(n):
        coeffs[i].append([])
    #Evaluating at theta = 0, which is the beginning of each of the pieces
    coeffs[i][0]=res['s'][i]
    try:
        coeffs[i][1]=res['v'][i]*((res['theta'][i+1]-res['theta'][i])/w)
        coeffs[i][2]=res['a'][i]/2*((res['theta'][i+1]-res['theta'][i])/w)**2
        coeffs[i][3]=res['j'][i]/6*((res['theta'][i+1]-res['theta'][i])/w)**3
    except:
        coeffs[i][1]=res['v'][i]*((res['theta'][0]-res['theta'][i]
            +2*np.pi)/w)
        coeffs[i][2]=res['a'][i]/2*((res['theta'][0]-res['theta'][i]
            +2*np.pi)/w)**2
        coeffs[i][3]=res['j'][i]/6*((res['theta'][0]-res['theta'][i]
            +2*np.pi)/w)**3
    #Time for some linear algebra
    #First construct the vector of all the knowns
    try:
        knowns=np.array([res['j'][i+1]*((res['theta'][i+1]
            -res['theta'][i])/w)**3-6*coeffs[i][3],
            res['a'][i+1]*((res['theta'][i+1]-res['theta'][i])/w)**2
            -(2*coeffs[i][2]+6*coeffs[i][3]),
            res['v'][i+1]*((res['theta'][i+1]-res['theta'][i])/w)
            -(coeffs[i][1]+2*coeffs[i][2]+3*coeffs[i][3]),
            res['s'][i+1]-(coeffs[i][0]+coeffs[i][1]+coeffs[i][2]
                +coeffs[i][3])])
    except:
        knowns=np.array([res['j'][0]*((res['theta'][0]
            -res['theta'][i]+2*np.pi)/w)**3-6*coeffs[i][3],
            res['a'][0]*((res['theta'][0]-res['theta'][i]+2*np.pi)/w)**2
            -(2*coeffs[i][2]+6*coeffs[i][3]),
            res['v'][0]*((res['theta'][0]-res['theta'][i]+2*np.pi)/w)
            -(coeffs[i][1]+2*coeffs[i][2]+3*coeffs[i][3]),
            res['s'][0]-(coeffs[i][0]+coeffs[i][1]+coeffs[i][2]
                +coeffs[i][3])])
    #Next resize for use with proper amount of coefficients
    knores=[]
    for j in range(8-n,8-4):
        knores.append(knowns[j])
    knores=np.array(knores)
    #Then multiply by the inverse of matrix found earlier
    try:
        found=np.dot(Minv,knores)
    except:
        found=Minv*knores[0][0]
    #add to the coefficient set
    for j in range(4,n):
        coeffs[i][j]=found[j-4]

#move results to database, reversing the order of coefficients to match polyval
header=[]
for i in range(n):
    header.append('C'+str(n-i))
results = pd.DataFrame(columns=header)
values=[]
for i in range(len(res)):
    values.append([])
    for j in range(n):
        values[i].append(coeffs[i][n-j-1])
    results = pd.DataFrame(np.insert(results.values, i, values[i], axis=0))
#setting meaningful header names
if n == 5:
    results = results.rename(columns={0:'C4',1:'C3',2:'C2',3:'C1',4:'C0'})
elif n == 6:
    results = results.rename(columns={0:'C5',1:'C4',2:'C3',3:'C2',4:'C1',
        5:'C0'})
elif n == 7:
    results = results.rename(columns={0:'C6',1:'C5',2:'C4',3:'C3',4:'C2',
        5:'C1',6:'C0'})
else:
    results = results.rename(columns={0:'C7',1:'C6',2:'C5',3:'C4',4:'C3',
        5:'C2',6:'C1',7:'C0'})

#new name for file to store results
filename=filename[:-4]+'cof.csv'

#Finally export the file
results.to_csv(filename, index=False)
