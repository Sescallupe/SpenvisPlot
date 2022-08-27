import numpy  as np
x=0.001*np.power(10,np.arange(21.)/5.)
y=np.exp(-x)
yint=y-np.exp(-10)
yint1=0.*yint
yint1[1:]=yint[0:-1]-yint[1:]
np.savetxt("user_table.txt", np.transpose(np.array([x,yint1])), "%.3e")
np.savetxt("arb_table.txt", np.transpose(np.array([x,y])), "%.3e")
