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
        jcof.iloc[i,j]*=((len(cof.iloc[i])-j-1)*(len(cof.iloc[i])-j-3)
                *(len(cof.iloc[i])-j-2))

#s formula
def smarch(i):
    i-=int(i/(2*np.pi))*2*np.pi
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
    i-=int(i/(2*np.pi))*2*np.pi
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
    i-=int(i/(2*np.pi))*2*np.pi
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
    i-=int(i/(2*np.pi))*2*np.pi
    if i>=res.loc[0,'theta'] and i<res.loc[1,'theta']:
        jerk=np.polyval(jcof.iloc[0],
                (i/(res.loc[1,'theta']-res.loc[0,'theta'])))
        jerk*=(w/(res.loc[1,'theta']-res.loc[0,'theta']))**3
    elif i<res.loc[2,'theta']:
        jerk=np.polyval(jcof.iloc[1],((i-res.loc[1,'theta'])
            /(res.loc[2,'theta']-res.loc[1,'theta'])))
        jerk*=(w/(res.loc[2,'theta']-res.loc[1,'theta']))**3
    elif i<res.loc[3,'theta']:
        jerk=np.polyval(jcof.iloc[2],((i-res.loc[2,'theta'])
            /(res.loc[3,'theta']-res.loc[2,'theta'])))
        jerk*=(w/(res.loc[3,'theta']-res.loc[2,'theta']))**3
    else:
        jerk=np.polyval(jcof.iloc[3],((i-res.loc[3,'theta'])
            /(res.loc[0,'theta']-res.loc[3,'theta']+2*np.pi)))
        jerk*=(w/(res.loc[0,'theta']-res.loc[3,'theta']+2*np.pi))**3
    return jerk

#Length of follower
l=1

#Prime Diameter
prime = 0.375

#Create a figure that this will be animated on
fig, ax = plt.subplots()
cam, = ax.plot([],[],'k-')
fol, = ax.plot([],[],'r-')
#fig, ax = plt.subplots(4,1)
#dis, = ax[0].plot([],[],'r.')
#vel, = ax[1].plot([],[],'r.')
#acc, = ax[2].plot([],[],'r.')
#jer, = ax[3].plot([],[],'r.')

#Full theta
thf=np.linspace(0,2*np.pi,256)

#Evaluating length for boundaries
standin=[0]*256
for j in range(256):
    standin[j]=smarch(thf[j])
raiser=min(standin)

#Cam init
def cinit():
    ax.set_xlim(-(max(standin)+raiser+prime),max(standin)+raiser+prime)
    ax.set_ylim(-(prime+raiser+max(standin)),prime+l+raiser+max(standin))
    ax.set_aspect('equal')
    return fol, cam,

#Cam rotation
x, y = [0]*256,[0]*256
xfol,yfol = [0]*4,[0]*4
def camrotate(i):
    for j in range(256):
        thi=thf[j]+i
        x[j]=(prime/2-raiser+smarch(thi))*np.sin(thf[j])
        y[j]=(prime/2-raiser+smarch(thi))*np.cos(thf[j])
    yfol=[prime/2+smarch(i),prime/2+smarch(i)+l,prime/2+smarch(i)+l,prime/2+smarch(i)+l]
    xfol=[0,0,1,-1]
    cam.set_data(x,y)
    fol.set_data(xfol,yfol)
    return fol, cam,

#Make still graphs
if False:
    seval=np.linspace(0,1,256)
    veval=np.linspace(0,1,256)
    aeval=np.linspace(0,1,256)
    jeval=np.linspace(0,1,256)
    for j in range(256):
        seval[j]=smarch(thf[j])
        veval[j]=vmarch(thf[j])
        aeval[j]=amarch(thf[j])
        jeval[j]=jmarch(thf[j])
    ax[0].plot(thf,seval,'k-')
    ax[1].plot(thf,veval,'k-')
    ax[2].plot(thf,aeval,'k-')
    ax[3].plot(thf,jeval,'k-')

#init functions
def sinit():
    ax[0].set_xlim(0,2*np.pi)
    return dis,
def vinit():
    ax[1].set_xlim(0,2*np.pi)
    return vel,
def ainit():
    ax[2].set_xlim(0,2*np.pi)
    return acc,
def jinit():
    ax[3].set_xlim(0,2*np.pi)
    return jer,

#Following along on the graphs
xdis,ydis=deque(maxlen=1),deque(maxlen=1)
def sfollow(i):
    ydis.appendleft(smarch(i))
    xdis.appendleft(i)
    dis.set_data(xdis,ydis)
    return dis,
xvel,yvel=deque(maxlen=1),deque(maxlen=1)
def vfollow(i):
    yvel.appendleft(vmarch(i))
    xvel.appendleft(i)
    vel.set_data(xvel,yvel)
    return vel,
xacc,yacc=deque(maxlen=1),deque(maxlen=1)
def afollow(i):
    yacc.appendleft(amarch(i))
    xacc.appendleft(i)
    acc.set_data(xacc,yacc)
    return acc,
xjer,yjer=deque(maxlen=1),deque(maxlen=1)
def jfollow(i):
    yjer.appendleft(jmarch(i))
    xjer.appendleft(i)
    vel.set_data(xjer,yjer)
    return jer,

def init():
    for j in range(4):
        ax[j].set_xlim(0,2*np.pi)
    return dis, vel, acc, jer

def march(i):
    ydis.appendleft(smarch(i))
    xdis.appendleft(i)
    dis.set_data(xdis,ydis)
    yvel.appendleft(vmarch(i))
    xvel.appendleft(i)
    vel.set_data(xvel,yvel)
    yacc.appendleft(amarch(i))
    xacc.appendleft(i)
    acc.set_data(xacc,yacc)
    yjer.appendleft(jmarch(i))
    xjer.appendleft(i)
    vel.set_data(xjer,yjer)
    return dis, vel, acc, jer,

gran=256
delay=w/(2*np.pi)*gran
#The animation
caman=ani.FuncAnimation(fig,camrotate,frames=np.linspace(0,2*np.pi,gran),
        interval = delay,init_func=cinit, blit=True)
#graphs=ani.FuncAnimation(fig,march,frames=np.linspace(0,2*np.pi,256),
#        init_func=init, blit=True)
#san=ani.FuncAnimation(fig,sfollow,frames=np.linspace(0,2*np.pi,256),
#        init_func=sinit, blit=True)
#van=ani.FuncAnimation(fig,vfollow,frames=np.linspace(0,2*np.pi,256),
#        init_func=vinit, blit=True)
#aan=ani.FuncAnimation(fig,afollow,frames=np.linspace(0,2*np.pi,256),
#        init_func=ainit, blit=True)
#jan=ani.FuncAnimation(fig,jfollow,frames=np.linspace(0,2*np.pi,256),
#        init_func=jinit, blit=True)
plt.show()

filename = filename[:-7]+'cam.mp4'
caman.save(filename)
