import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.integrate as integrate


def find_fourier_coef(a,phi,num=100,bc_type="dirichlet",bc_vals=(0,0)):
    """Generates num coefficients of the Fourier series of phi for the
    given interval length and boundary conditions. """
    if bc_type == "dirichlet":
        b = np.zeros(num)
        for n in range(1,num):
            func = lambda x:phi(x)*np.sin(n*x*np.pi/a)
            b[n] = (2/a)*integrate.quad(func,0,a)[0]

    ##TODO test this
    ##TODO make work for Robin/ different bcvals
    if bc_type == "neumann":
        b = np.zeros(num)
        for n in range(num):
            func = lambda x:phi(x)*np.cos(n*x*np.pi/a)
            b[n] = (2/a)*integrate.quad(func,0,a)[0]

    return b

def plot_diffusion(D,a,phi,t_min=0,t_max=1,bc_type="dirichlet"):
    """Uses Fourier Series to plot """
    ##creating x and t domains
    num_points = 100
    x = np.linspace(0,a,num_points)
    t = np.linspace(t_min,t_max,num_points)
    X,T = np.meshgrid(x,t)

    ##finding fourier coefficients
    num_sin =100
    b = find_fourier_coef(a,phi,num = num_sin,bc_type="dirichlet")
    
    ##creating data
    data =np.zeros(X.shape)
    for n in range(num_sin):
        data += b[n]*np.sin(n*np.pi*X/a)*np.exp(-D*((n*np.pi)**2)*T)

    
    ##initial plot 
    fig,ax = plt.subplots()
    ax.set(xlim=(0,a),ylim =(0,2))
    line = ax.plot(x,data[0,:])[0]
    
    def update(i):
        line.set_ydata(data[i,:])
    anim = animation.FuncAnimation(fig,update,interval =40,frames=len(t)-1)
    return anim

def plot_diffusion_alpha(alpha,D,a,phi,t_min=0,t_max=1,bc_xmin=0,bc_xmax=0):
    """Uses Fourier Series to plot """
    ##creating x and t domains
    num_points = 100
    x = np.linspace(0,a,num_points)
    t = np.linspace(t_min,t_max,num_points)
    X,T = np.meshgrid(x,t)

    ##finding fourier coefficients
    num_sins = 100
    b = np.zeros(num_sins)
    for n in range(1,num_sins):
        func = lambda x:phi(x)*np.sin(n*x*np.pi/a)
        b[n] = (2/a)*integrate.quad(func,0,a)[0]
    
    ##creating data
    data =np.zeros(X.shape)
    for n in range(num_sins):
        data += b[n]*np.sin(n*np.pi*X/a)*np.exp((alpha-(D*(n*np.pi)**2)/a**2)*T)
    print("data generated")
    ##initial plot 
    fig,ax = plt.subplots()
    ax.set(xlim=(0,a),ylim =(0,2))
    line = ax.plot(x,data[0,:])[0]
    
    def update(i):
        line.set_ydata(data[i,:])
        
    anim = animation.FuncAnimation(fig,update,interval =40,frames=len(t)-1)
    return anim

def plot_inhomogenous_diffusion(D,a,phi,f,t_min=0,t_max=1,bc_type="dirichlet"):
    """Uses Fourier Series to plot """
    ##creating x and t domains
    num_points = 100
    x = np.linspace(0,a,num_points)
    t = np.linspace(t_min,t_max,num_points)
    X,T = np.meshgrid(x,t)

    ##finding fourier coefficients
    num_sin =100
    phi_series = find_fourier_coef(a,phi,num = num_sin,bc_type="dirichlet")
    f_series = find_fourier_coef(a,f,num=num_sin,bc_type="dirichlet")

    ##creating data
    ##TODO need to work out integrals from 0 to t of f_n functions, as f is time dependant 
    data =np.zeros(X.shape)
    for n in range(num_sin):
        gamma_n = D*((n*np.pi)**2)
        u_n = phi_series[n]+integrate.quad(f_n*np.exp(gamma_n),0,)
        data += u_n*np.sin(n*np.pi*X/a)*np.exp(-gamma_n*T)

    
    ##initial plot 
    fig,ax = plt.subplots()
    ax.set(xlim=(0,a),ylim =(0,2))
    line = ax.plot(x,data[0,:])[0]
    
    def update(i):
        line.set_ydata(data[i,:])
    anim = animation.FuncAnimation(fig,update,interval =40,frames=len(t)-1)
    return anim

# def dirichlet_solution(n,X,T,a,D):
#     return np.sin(n*np.pi*X/a)*np.exp(-D*((n*np.pi)**2)*T)
# def neumann_solution(n,X,T,a,D):
#     pass
class PDE_Toolkit():
    """Finds fourier decompositions and time dependant evolutions
    of given initial condition and equation type. Currently supported
    equation types """
    # SOL_DICT = {
    #     "dirichlet":dirichlet_solution,
    #     "neumann": neumann_solution
    #     }
    def __init__(self,initial_condition,equation="diffusion",bc_type =["d",[0,0]]):
        self.ic = initial_condition
        self.eq = "diffusion"

    def find_fourier_coeff(self):
        pass

    def generate_data(self,t_min,t_max,dt=0.05,save_rate=1):
        T = t_max-t_min
        num_time_steps = int(T//dt)
        num_space_steps = 100
        self.data = np.zeros()
    
    def animate_data(self,t_min,t_max,dt=0.05,save_rate=1):
        fig,ax = plt.subplots()
        ax.set(xlim=(0,a),ylim =(0,2))
        line = ax.plot(x,data[0,:])[0]
        
        def update(i):
            line.set_ydata(data[i,:])
        anim = animation.FuncAnimation(fig,update,interval =40,frames=len(t)-1)
        return anim