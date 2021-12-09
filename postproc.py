import numpy as np
import matplotlib.pyplot as plt
import matplotlib  

matplotlib.rcParams.update({'font.size': 14}) #Default Font Size
matplotlib.rcParams['text.usetex'] = True

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def sin(a):
    return np.sin(a)

def cos(a):
    return np.cos(a)

def get_data(fname,col,nx,ny):
    q_data = np.genfromtxt(fname, skip_header=3)
    x1 = q_data[:,0]
    y1 = q_data[:,1]
    q1 = q_data[:,col-1]
    xq = np.reshape(x1,(nx,ny),order='F') # order F for fortran is important for proper reshape
    yq = np.reshape(y1,(nx,ny),order='F')
    data  = np.reshape(q1,(nx,ny),order='F')
    return xq,yq,data

def plot_contour(xq,yq,data,title,fname):
    plt.figure(figsize=(8,8))
    plt.contourf(xq,yq,data,cmap='YlOrRd')
    plt.colorbar()
    plt.title(title)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.savefig(fname)
    plt.show()  # show the plot     
    return

def plot_data(a,b,title,fname):
    plt.figure(figsize=(8,8))
    plt.plot(a,b)
    plt.title(title)
    plt.xlabel(r'Residual')
    plt.ylabel(r'Iterations')
    plt.savefig(fname)
    plt.show()  # show the plot     
    return

xq,yq,data = get_data('data.dat',3,9,9)
plot_contour(xq,yq,data,'Contour','contour.png')

logdata = np.genfromtxt('log.dat', skip_header=2)
plot_data(logdata[:,0],logdata[:,1],'Residual vs Iterations','residual.png')
