# Theory

Biofilm.jl simulates a one-dimensional biofilm within a stirred tank reactor.  The dependent variables include the 
- tank particulates (biomass) concentration(s) ``X``,
- tank solute concentrations ``S``,
- biofilm particulate (biomass) volume fractions ``P_b``,
- biofilm solute concentrations ``S_b``, and
- biofilm thickness ``L_f``. 

## Tank Equations
### Particulates
The governing equation describing the particulate concentrations in the tank environment is
```math
\frac{d X_{t,j}}{dt} = \mu_j(\mathbf{S}_t) X_{t,j} - \frac{Q X_{t,j}}{V} + \frac{v_\mathrm{det} A X_{b,j}(L_f)}{V} + \mathrm{src}_{X,j}
```
for ``j=1,\dots,N_x``, where ``t`` is time, ``\mu_j(\mathbf{S}_t)`` is the growthrate of the ``j^\mathrm{th}`` particulate, ``Q`` is the flowrate, ``V`` is the volume of the tank, ``v_\mathrm{det}=K_\mathrm{det} L_f^2`` is the detachment velocity, ``A`` is the area of the biofilm, ``X_{b,j}(L_f)`` is the ``j^\mathrm{th}`` particulate concentration at the top of the biofilm, and ``\mathrm{src}_{X,j}`` is the source term for the ``j^\mathrm{th}`` particulate. 

The terms on the right-hand-side (RHS) are 
1) the growth of the particulate in the tank, 
2) transport due to flow out of the tank, 
3) transfer of particulates from the biofilm to the tank due to detachment, and
4) source term.

### Solutes
The governing equation describing the solute concentrations in the tank environment is
```math
\frac{d S_{t,k}}{dt} = -\sum_{j=1}^{N_x} \frac{\mu_j(\mathbf{S}_t) X_{t,j}}{Y_{j,k}} + \frac{Q S_{\mathrm{in},k}}{V} - \frac{Q S_{t,k}}{V} + \frac{A S_{\mathrm{flux},k}}{V} + \mathrm{src}_{S,k}
```
for ``k=1,\dots,N_s``, where ``S_{\mathrm{flux},k}`` is the flux of solutes from the biofilm into the tank, and ``\mathrm{src}_{S,k}`` is the source term for the ``k^\mathrm{th}`` solute. 

The terms on the right-hand-side (RHS) are 
1) consumption of solutes due to the growth of the particulate in the tank, 
2) transport due to flow into the tank, 
3) transport due to flow out of the tank,
4) transfer of solutes into the biofilm due to diffusion, and
5) source term.
   
## Biofilm Equations
### Particulates
The governing equations describing the biofilm environment are
```math
\frac{d P_{b,j,i}}{dt} = 
\mu_j(\mathbf{S}_{b,i}) P_{b,j,i} 
- \frac{d v_i P_{b,j,i}}{dz} 
+ \frac{\mathrm{src}_{X,j,i}}{\rho_j}
```
for ``j=1,\dots,N_x`` and ``i=1,\dots,N_z``. Where ``P_{b,j,i}`` is the ``j^\mathrm{th}`` particulate at the ``i^\mathrm{th}`` grid point within the biofilm. 

The terms on the right-hand-side (RHS) are 
1) the growth of the particulate in the biofilm, 
2) transport through the biofilm due to the growth velocity ``v_i``, and 
3) source term of particulate at ``i^\mathrm{th}`` location in biofilm.

The growth velocity ``v_i`` is the rate of flow through the biofilm due to growth deeper within the biofilm and is defined with
```math
v_i=  \int_{z=0}^{z_i}{\sum_{j=1}^{N_x} \frac{1}{P_\mathrm{tot}}\left(\mu_j(\mathbf{S}_{b,i}) P_{b,j,i} + \frac{\mathrm{src}_{X,j}}{\rho_j}\right) ~dz}
```
where ``P_\mathrm{tot}=\sum_{j=1}^{N_x}{P_{b,j}}``
   
### Solutes
```math
\frac{d S_{b,k,i}}{dt} = 
D_{e,k}\frac{d^2 S_{b,k,i}}{dz^2} 
- \sum_{j=1}^{N_x} \frac{\mu_j(\mathbf{S}_{b,i}) X_{b,j,i}}{Y_{j,k}}
+ \mathrm{src}_{S,k,i}
```
for ``k=1,\dots,N_s`` and ``i=1,\dots,N_z``.

The terms on the right-hand-side (RHS) are 
1) diffusion of solutes in the biofilm,
2) consumption of solutes due to the growth of the particulate in the biofilm, and
3) source term of solute at ``i^\mathrm{th}`` location in biofilm.

The diffusion term with a second derivative w.r.t. ``z`` requires boundary conditions at the top and bottom of the biofilm.  A zero-flux (zero first-derivative) condition is used at the bottom of the biofilm.  At the top of the biofilm the diffusion through the boundary layer is matched with the diffusion into the biofilm, i.e.,
```math
D_{aq,k}\frac{d^2 S_{k}}{dz^2} = D_{e,k}\frac{d^2 S_{b,k,N_z}}{dz^2} 
```
for ``k=1,\dots,N_s``, where ``D_e`` is the diffusion coefficient in the biofilm and ``D_\mathrm{aq}`` is the diffusion coefficient in the boundary layer. 

### Biofilm Thickness
The thickness of the biofilm ``L_f`` is described by 
```math
\frac{d L_f}{dt} = v_{N_z} - v_\mathrm{det}
```
where the first term on the RHS is the growth velocity at the top of the biofilm (see [Biofilm Particulates](#Particulates-2)) and the second term is the detachment velocity modeled with ``v_\mathrm{det}=K_\mathrm{det} L_f^2``