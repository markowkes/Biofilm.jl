# Case Parameters

Different biofilms can be modeled by modify the Case1.jl file.  This is a text file and can be opened and edited using Notepad (Windows), textEdit (Mac), or similar. 

The Case file defines a param struct with the following fields

### Simulation Parameters
- **Title** - Description of the case, used, e.g., on the title of plots 
- **tFinal** - Simulation is performed from ``t=0`` to ``t=t_\mathrm{Final}`` [days]
- **tol** - Tolerance used for differential equation solver.  Solution will have error ``\le`` tol.
- **outPeriod** - Period in days between outputs during simulation.  Note that the solver will take smaller timesteps then outPeriod and the entire solution will be available when the solver completes (that is outPeriod has no impact on the solution). [days]
- **plotPeriod** - *[default=outPeriod]* - Period in days between plot renderings during the simulation.  For long simulations you can increase plotPeriod to reduce the number of plots that are created and speedup the simulation. Note: plotPeriod is required to be a multiple of outPeriod. [days]

### Particulate Parameters
- **XNames** - Names of the particulates used in plots and other outputs
- **Xto** - Initial condition of particulate concentration(s) in tank environment [g/m``^3``]
- **Pbo** - Initial condition of particulate volume fraction(s) in biofilm [-]
- **rho** - Particulate densities [g/m``^3``]
- **Kdet** - Particulate detachment coefficient [1/(m days)]
- **src** - Array of source terms for each particulate.  Each source term can be a function of solute or particulate concentrations.  See Case 4 for an example of using src. [g/(m``^3`` s)]
- **mu** - Growthrates for each particulate.  Each growthrate can be a function of solute and particulate concentrations, biofilm thickness, time, position within biofilm, and parameters.  See example for additional details. [1/days]
- **discontinuityPeriod** - [Optional] - Period in days between discontinuities in mu or src.  If mu or src are continuous do not set this optional parameter. [days]

### Solute Parameters
- **SNames** - Names of the solutes used in plots and other outputs
- **Sto** - Initial condition of solute concentration(s) in tank environment [g/m``^3``]
- **Sbo** - Initial condition of solute concentrations(s) in biofilm [g/m``^3``]
- **Sin** - Vector of functions that provide the inflow solute concentrations [g/m``^3``].
- **Yxs** - Array of biomass yield coefficients [g``_\mathrm{X}``/g``_\mathrm{S}``]
  ```math
  Y_{xs} = \begin{bmatrix}
  \frac{\Delta X_1}{\Delta S_1} & \frac{\Delta X_1}{\Delta S_2} & \cdots & \frac{\Delta X_1}{\Delta S_N}\\[0.5em]
  \frac{\Delta X_2}{\Delta S_1} & \frac{\Delta X_2}{\Delta S_2} & \cdots & \frac{\Delta X_2}{\Delta S_N}\\[0.5em]
  \vdots & \vdots & \ddots & \vdots \\[0.5em]
  \frac{\Delta X_M}{\Delta S_1} & \frac{\Delta X_M}{\Delta S_2} & \cdots & \frac{\Delta X_M}{\Delta S_N}\\
  \end{bmatrix}
  ```
  for ``N`` solutes and ``M`` particulates.  If ``X_i`` does not depend on ``S_j``, set ``\displaystyle \frac{\Delta X_i}{\Delta S_j}=0``.
- **Dt** - Aqueous diffusion coefficient of solutes through a boundary layer within the tank at the top of the biofilm.  If you do not want a boundary layer set **LL=0**, then these values will not impact the solution. [m``^2``/day]
- **Db** - Effective diffusion coefficient of solutes through a biofilm [m``^2``/day]

### Tank Parameters
- **V** - Volume of tank [m``^3``]
- **S** - Surface area of biofilm [m``^2``]
- **Q** - Flowrate through tank [m``^3``/day]

### Biofilm Parameters
- **Nz** - Number of grid points used to discretize the biofilm
- **Lfo** - Initial condition for biofilm thickness [m]
- **LL** - Thickness of boundary layer at top of biofilm [m]


