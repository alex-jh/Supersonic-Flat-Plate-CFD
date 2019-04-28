## Hello world

a $k_{n+1}$

$$k_{n+1}$$

I have created a class that contains all the data specific to this problem, so later I can reuse code for other configurations. So let's call the class SupersonicFlatPlate.

The first step is to initialize the relevant parameters:

$$M_\infty = 4$$
$$Plate Length = 0.00001m$$
$$a_\infty = 340.28m/s$$
$$P_\infty = 10132503N/m^2$$
$$T_\infty = 288.16K$$
$$T_{wall}/T_\infty = 1$$
$$\gamma = 0.4$$
$$P_r = 0.71$$
$$\mu_0 = 1.7894\cdot 10^{-5}kg/{(m\cdot s)}$$
$$T_0 = 288.16K$$
$$R = 287J/{(Kg\cdot K)}$$

The plate length is so small to make the solution converge faster, and to reduce the size of the mesh.

## Calculation of Step Sizes in Space

The first thing we need to do is define the grid size in the streamwise direction. I´m going to use a grid of 100 x 100. With this we can calculate the step size in the x direction as:

$$IMAX=100$$
$$\Delta  x = \frac{LengthPlate}{IMAX}$$

For the y direction we also have 100 grids, but we need to guess the size of the domain. For this, we can assume it is going to be five times the size of the boundary layer, as predicted by the Blasius equation, using the Reynolds number.

$$L_{vert} = 5 \cdot \delta$$
$$\delta = \frac{5 \cdot LengthPlate}{\sqrt{R_e}}$$
$$R_e = \rho u LengthPlate / \mu$$

from which we can calculate the step size in the y drection as:

$$\Delta  x = \frac{L_{vert}}{JMAX}$$

We could verify these step size by calculating the Reynolds numbers and set some appropiate values. The requirements in the y direction should be stricter than in the y direction, because it has stronger gradients.

## Step Size in Time
To calculate the step size, we are going to use the Courant-Friedrichs-Lewy (CFL) criterion, which it can be defined as follows:

$$(\Delta t_{CFL})_{i,j} = [\frac{|u_{i,j}}{\Delta x} + \frac{|v_{i,j}}{\Delta y} + a_{i,j}\sqrt{\frac{1}{\Delta x^2} + \frac{1}{\Delta y^2}} + 2v'_{i,j}(\frac{1}{\Delta x^2} + \frac{1}{\Delta y^2})]^{-1}$$

where

$$v'_{i,j} = max[\frac{4/3\mu_{i,j}(\gamma \mu_{i,j}/P_r)}{\rho_{i,j}}]$$

Then

$$\Delta t = min[K(\Delta t_{CFL})_{i,j}]$$

and K is the Courant factor that acts as a "fudge factor" to make sure the solutin remains stable.

### Viscosity
To calculate the viscosity, assuming perfect gas, we are going to use Sutherland's law.

$$\mu=\mu_0*(\frac{T}{T_o})^{3/2}\frac{T_0+110}{T+110}$$

## Initial Condition

We can set the intial condition equal to the freestream flow except over the plate, where we set the the wall Temperature and a no-slip condition

$$u = v = 0.0$$

The pressure will be calculated by extrapolating the pressure near the edge:

$$P(i,0) = 2P(i, 1) - P(i, 2)$$

## Boundary Conditions

Having specified the inital conditions, we need to enforce boundary conditions.

### (0,0)
In the leading edge, which is the point (0, 0), we enforce the no-slip condition, and the temperatur and pressure are asumed to take the freestream conditions.

### Left and Upper Boundary
In the left and upper boundary, the flow is assumed to take the freestream conditions, and the y component of the velocity is assumed to be equal to zero.

### Surface of the Plate
In the surface of the plate, we add a the no-slip condition, a temperature equal to the wall temperature and we extrapolate the pressure, in the same way as we have done for the initial conditions.

$$u = v = 0.0$$
$$T = T_{wall}$$
$$P(i,0) = 2P(i, 1) - P(i, 2)$$

### Right Side
All properties in the right side are calculated by extrapolation of interior points. For example, the velocity u is:

$$u(imax, j) = 2u(imax - 1, j) - u(imax - 2, j)$$

## MacCormack Algorithm

~~~~
$M_\inf$=4
~~~~
