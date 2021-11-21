# PDE Toolkit 
## Overview
The PDE Toolkit class allows you to quickly generate numerical fourier series, sample data for a given xrange and trange and animate this data for the wave and diffusion equations with Neumann and Dirichlet boundary conditions. 
## The Maths 
PDE Toolkit supports 4 equation types in 1D:
### Homogeneous Diffusion Equation
$$u_t(x,t)=Du_{xx}(x,t)$$
This has general solution
$\begin{align}
x=1
\end{align}$
### Inhomogeneous Diffusion Equation
$$u_t(x,t)=Du_{xx}(x,t)+f(x,t)$$
### Homogeneous Wave Equation
$$u_{tt}(x,t)=c^2u_{xx}(x,t)$$
### Homogeneous Diffusion Equation
$$u_{tt}(x,t)=c^2u_{xx}(x,t)+f(x,t)$$

