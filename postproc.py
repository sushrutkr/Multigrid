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

def plot_data(a,b,title,fname,label):
    plt.plot(a,b,label=label)
    plt.title(title)
    plt.legend()
    plt.ylabel(r'Residual')
    plt.xlabel(r'Iterations')

    # plt.xscale("log")
    plt.yscale("log")
    
    return

xq,yq,data = get_data('../Code/data.dat',3,65,65)
plot_contour(xq,yq,data,'Contour','contour.png')

logdata_mg = np.genfromtxt('../Code/log_mg.dat', skip_header=2)
logdata_gs = np.genfromtxt('../Code/log_gs.dat', skip_header=2)

plt.figure(1)
plot_data(logdata_mg[:,0],logdata_mg[:,1],'Residual vs Iterations','residual.png', 'Multigrid')
plot_data(logdata_gs[:,0],logdata_gs[:,1],'Residual vs Iterations','residual.png', 'Gauss-Seidel')

plt.savefig('residual.eps')
plt.show()  # show the plot 