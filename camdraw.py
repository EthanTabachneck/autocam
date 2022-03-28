import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from collections import deque
import pandas as pd

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

filename = filename[:-4]+'cof.csv'
coeff=pd.read_csv(filename)

#Convert datafram to float
for i in range(len(res)):
    res.loc[i,'w']=float(res.loc[i,'w'])
    res.loc[i,'theta'] = float(res.loc[i,'theta'])
    res.loc[i,'s'] = float(res.loc[i,'s'])
    res.loc[i,'v'] = float(res.loc[i,'v'])
    res.loc[i,'a'] = float(res.loc[i,'a'])
    res.loc[i,'j'] = float(res.loc[i,'j'])
    try:
        coeff.loc[i,'C7']=float(coeff.loc[i,'C7'])
    except:
        a='not 8 coefficients'
    try:
        coeff.loc[i,'C6']=float(coeff.loc[i,'C6'])
    except:
        a='not 7 coefficients'
    try:
        coeff.loc[i,'C5']=float(coeff.loc[i,'C5'])
    except:
        a='not 6 coefficients'
    coeff.loc[i,'C4']=float(coeff.loc[i,'C4'])
    coeff.loc[i,'C3']=float(coeff.loc[i,'C3'])
    coeff.loc[i,'C2']=float(coeff.loc[i,'C2'])
    coeff.loc[i,'C1']=float(coeff.loc[i,'C1'])
    coeff.loc[i,'C0']=float(coeff.loc[i,'C0'])

#conversion to radians
for i in range(len(res)):
    res.loc[i,'w']=np.radians(res.loc[i,'w'])
    res.loc[i,'theta']=np.radians(res.loc[i,'theta'])

#Setting w to w
w=res['w'][0]

#making figures
fig1, axs = plt.subplots(4,1,constrained_layout=True)

#Gathering coefficients
cof=coeff.iloc[:,:].copy()
vcof=coeff.iloc[:,:-1].copy()
acof=coeff.iloc[:,:-2].copy()
jcof=coeff.iloc[:,:-3].copy()
for i in range(len(cof)):
    for j in range(len(vcof.iloc[0])):
        vcof.iloc[i,j]*=(len(cof.iloc[i])-j-1)
    for j in range(len(acof.iloc[0])):
        acof.iloc[i,j]*=(len(cof.iloc[i])-j-1)*(len(cof.iloc[i])-j-2)
    for j in range(len(jcof.iloc[0])):
        jcof.iloc[i,j]*=(len(cof.iloc[i])-j-1)*(len(cof.iloc[i])-j-3)*(len(cof.iloc[i])-j-2)

#s formula
def smarch(i):
    if i>=res.loc[0,'theta'] and i<res.loc[1,'theta']:
        s=np.polyval(cof.iloc[0],((i-res.loc[0,'theta'])
            /(res.loc[1,'theta']-res.loc[0,'theta'])))
    elif i<res.loc[2,'theta']:
        s=np.polyval(cof.iloc[1],((i-res.loc[1,'theta'])
            /(res.loc[2,'theta']-res.loc[1,'theta'])))
    elif i<res.loc[3,'theta']:
        s=np.polyval(cof.iloc[2],((i-res.loc[2,'theta'])
            /(res.loc[3,'theta']-res.loc[2,'theta'])))
    else:
        s=np.polyval(cof.iloc[3],((i-res.loc[3,'theta'])
            /(res.loc[0,'theta']-res.loc[3,'theta']+2*np.pi)))
    return s

#v formula
def vmarch(i):
    if i>=res.loc[0,'theta'] and i<res.loc[1,'theta']:
        v=np.polyval(vcof.iloc[0],(i/(res.loc[1,'theta']-res.loc[0,'theta'])))
        v*=(w/(res.loc[1,'theta']-res.loc[0,'theta']))
    elif i<res.loc[2,'theta']:
        v=np.polyval(vcof.iloc[1],((i-res.loc[1,'theta'])
            /(res.loc[2,'theta']-res.loc[1,'theta'])))
        v*=(w/(res.loc[2,'theta']-res.loc[1,'theta']))
    elif i<res.loc[3,'theta']:
        v=np.polyval(vcof.iloc[2],((i-res.loc[2,'theta'])
            /(res.loc[3,'theta']-res.loc[2,'theta'])))
        v*=(w/(res.loc[3,'theta']-res.loc[2,'theta']))
    else:
        v=np.polyval(vcof.iloc[3],((i-res.loc[3,'theta'])
            /(res.loc[0,'theta']-res.loc[3,'theta']+2*np.pi)))
        v*=(w/(res.loc[0,'theta']-res.loc[3,'theta']+2*np.pi))
    return v

#a formula
def amarch(i):
    if i>=res.loc[0,'theta'] and i<res.loc[1,'theta']:
        a=np.polyval(acof.iloc[0],(i/(res.loc[1,'theta']-res.loc[0,'theta'])))
        a*=(w/(res.loc[1,'theta']-res.loc[0,'theta']))**2
    elif i<res.loc[2,'theta']:
        a=np.polyval(acof.iloc[1],((i-res.loc[1,'theta'])
            /(res.loc[2,'theta']-res.loc[1,'theta'])))
        a*=(w/(res.loc[2,'theta']-res.loc[1,'theta']))**2
    elif i<res.loc[3,'theta']:
        a=np.polyval(acof.iloc[2],((i-res.loc[2,'theta'])
            /(res.loc[3,'theta']-res.loc[2,'theta'])))
        a*=(w/(res.loc[3,'theta']-res.loc[2,'theta']))**2
    else:
        a=np.polyval(acof.iloc[3],((i-res.loc[3,'theta'])
            /(res.loc[0,'theta']-res.loc[3,'theta']+2*np.pi)))
        a*=(w/(res.loc[0,'theta']-res.loc[3,'theta']+2*np.pi))**2
    return a

#j formula
def jmarch(i):
    if i>=res.loc[0,'theta'] and i<res.loc[1,'theta']:
        j=np.polyval(jcof.iloc[0],(i/(res.loc[1,'theta']-res.loc[0,'theta'])))
        j*=(w/(res.loc[1,'theta']-res.loc[0,'theta']))**3
    elif i<res.loc[2,'theta']:
        j=np.polyval(jcof.iloc[1],((i-res.loc[1,'theta'])
            /(res.loc[2,'theta']-res.loc[1,'theta'])))
        j*=(w/(res.loc[2,'theta']-res.loc[1,'theta']))**3
    elif i<res.loc[3,'theta']:
        j=np.polyval(jcof.iloc[2],((i-res.loc[2,'theta'])
            /(res.loc[3,'theta']-res.loc[2,'theta'])))
        j*=(w/(res.loc[3,'theta']-res.loc[2,'theta']))**3
    else:
        j=np.polyval(jcof.iloc[3],((i-res.loc[3,'theta'])
            /(res.loc[0,'theta']-res.loc[3,'theta']+2*np.pi)))
        j*=(w/(res.loc[0,'theta']-res.loc[3,'theta']+2*np.pi))**3
    return j

#making sure it is not plotting past 2pi
axs[0].set_xlim(0,2*np.pi)
axs[1].set_xlim(0,2*np.pi)
axs[2].set_xlim(0,2*np.pi)
axs[3].set_xlim(0,2*np.pi)

#Full theta
thf=np.linspace(0,2*np.pi,256)

#plots
#Starting with making arrays of curve values
seval = np.linspace(0,1,256)
aeval = np.linspace(0,1,256)
veval = np.linspace(0,1,256)
jeval = np.linspace(0,1,256)
for k in range(256):
    seval[k] = smarch(thf[k])
    veval[k] = vmarch(thf[k])
    aeval[k] = amarch(thf[k])
    jeval[k] = jmarch(thf[k])
axs[0].plot(thf,seval,'k-')
axs[1].plot(thf,veval,'k-')
axs[2].plot(thf,aeval,'k-')
axs[3].plot(thf,jeval,'k-')

prime = 0.375

#cam drawing
fig2, ax = plt.subplots()
x=np.linspace(0,1,256)
y=np.linspace(0,1,256)
standin=[0]*256
for j in range(256):
    standin[j]=smarch(thf[j])
raiser=min(standin)
for j in range(256):
    x[j]=(prime/2-raiser+smarch(thf[j]))*np.sin(thf[j])
    y[j]=(prime/2-raiser+smarch(thf[j]))*np.cos(thf[j])
ax.plot(x,y,'k-')
ax.set_aspect('equal')
plt.show()

if True:
    f1fil='images'+filename[4:-7]+'curves.png'
    f2fil='images'+filename[4:-7]+'cam.png'
    fig1.savefig(f1fil,bbox_inches='tight')
    fig2.savefig(f2fil,bbox_inches='tight')
