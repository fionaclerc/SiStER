% SiStER_Input_File_coulomb_cliff_GEOM
%
% sets up rectangular ice cliff exposed by linear thinning of 
% buttressing ice shelf. F. Clerc 2018

% reference model
% DURATION OF SIMULATION AND FREQUENCY OF OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%
dt_out = 1;% output files every "dt_out" iterations

D = H-h;
ha = .2e3;% height of sticky air layer
% DOMAIN SIZE AND GRIDDING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xsize=3e3;
ysize=H + ha;

% gridding- from 0 to GRID.x(1), grid size is GRID.dx(1)
% from GRID.x(1) to GRID.x(2), grid size is GRID.dx(1) etc...
% same for y

GRID.dx(1)=dxx;
GRID.x(1)=xsize;
GRID.dy(1) = dyy;
GRID.y(1) = H + ha;
nx = length([0:GRID.dx:GRID.x]);

% LAGRANGIAN MARKERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mquad=8; % number of markers in the smallest quadrant
Mquad_crit=4; % minimum number of markers allowed in smallest quadrant (for reseeding)

% GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nphase=5; % number of phases

% phase 1
% air
GEOM(1).type=1; % 1 = layer (then specify top and bot) or 2 = circle % 1 = layer (then specify top and bot) or 2 = circle (then specify center and radius)
GEOM(1).top=0;
GEOM(1).bot=H + ha + 0.5;

% phase 2
% ice cliff
GEOM(2).type=3; % 1 = layer (then specify top and bot) or 2 = circle % 1 = layer (then specify top and bot) or 2 = circle (then specify center and radius)
GEOM(2).top=ha;
GEOM(2).bot=H + ha;
GEOM(2).left=0e3;
GEOM(2).right=xsize/2 + 0.5;

% phase 3
% water
GEOM(3).type=3; % 1 = layer (then specify top and bot) or 2 = circle % 1 = layer (then specify top and bot) or 2 = circle (then specify center and radius)
GEOM(3).top=ha + h-0.5;
GEOM(3).bot=H + ha;
GEOM(3).left=xsize/2;
GEOM(3).right=xsize;

% phase 4
% air
GEOM(4).type=3; % 1 = layer (then specify top and bot) or 2 = circle % 1 = layer (then specify top and bot) or 2 = circle (then specify center and radius)
GEOM(4).top=ha;
GEOM(4).bot=ha + h;
GEOM(4).left=xsize/2;
GEOM(4).right=xsize;

% phase 5
% ice shelf
GEOM(5).type=3; % 1 = layer (then specify top and bot) or 2 = circle % 1 = layer (then specify top and bot) or 2 = circle (then specify center and radius)
GEOM(5).top=ha;
GEOM(5).bot=ha + H;
GEOM(5).left=xsize/2;
GEOM(5).right=xsize;

% MATERIAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% creep laws of the form: pre^(-1/n)*epsII^((1-n)/n)*exp(E/(nRT))
% harmonically averaging diffusion creep, dislocation creep 
% (and plastic creep to simulate brittle failure)
visc_air = 1e7;
visc_water = 1e7;
% phase 1 "sticky air"
MAT(1).phase=1;
% density parameters
MAT(1).rho0=0.01;
MAT(1).alpha=0;
% thermal parameters
MAT(1).k=3;
MAT(1).cp=1000;
% elasticity 
MAT(1).G=1e18;
% diffusion creep parameters
MAT(1).pre_diff=.5/visc_air;
MAT(1).Ediff=0;
%MAT(1).ndiff=n_flow;
MAT(1).ndiff = 1;
% dislocation creep parameters
MAT(1).pre_disc=.5/visc_air;
MAT(1).Edisc=0;
%MAT(1).ndisc=n_flow;
MAT(1).ndisc = 1;
% plasticity
MAT(1).mu=0;
MAT(1).mumin=0;
MAT(1).Cmax=0.01e12;
MAT(1).Cmin=0.01e12;
MAT(1).ecrit=0.1;

MAT_end(1) = MAT(1);


% phase 2 "ice"
MAT(2).phase=2;
% density parameters
MAT(2).rho0=900;%kg/m^3
MAT(2).alpha=0;
% thermal parameters
MAT(2).k=2.1;%W/(m*K)
MAT(2).cp=2093;% J/(K*kg)
% elasticity 
MAT(2).G=Gice;
% diffusion creep parameters
MAT(2).pre_diff=A;% 1/A should be more 1.16e-26
MAT(2).Ediff=0;
MAT(2).ndiff=n_flow;
% dislocation creep parameters
MAT(2).pre_disc=A;
MAT(2).Edisc=0;
MAT(2).ndisc=n_flow;
% plasticity
MAT(2).mu=0;
MAT(2).mumin=0;
MAT(2).Cmax=tau_yield;
MAT(2).Cmin=tau_yield;
MAT(2).ecrit=0.1;

MAT_end(2) = MAT(2);


% phase 3 "ice -> water"
%MAT(3) = MAT(2);
MAT(3).phase = 3;


MAT_end(3).phase=3;
% density parameters
MAT_end(3).rho0=1000;
MAT_end(3).alpha=0;
% thermal parameters
MAT_end(3).k=5;
MAT_end(3).cp=4182;% J/(K*kg)
% elasticity 
MAT_end(3).G=1e18;
% diffusion creep parameters
MAT_end(3).pre_diff=.5/visc_water;
MAT_end(3).Ediff=0;
%MAT_end(3).ndiff=n_flow;
MAT_end(3).ndiff=1;

% dislocation creep parameters
MAT_end(3).pre_disc=.5/visc_water;
MAT_end(3).Edisc=0;
%MAT_end(3).ndisc=n_flow;
MAT_end(3).ndisc=1;

% plasticity
MAT_end(3).mu=0;
MAT_end(3).mumin=0;
MAT_end(3).Cmax=400e9;
MAT_end(3).Cmin=400e9;
MAT_end(3).ecrit=0.1;
MAT(3) = MAT_end(3);
%
% phase 4 "ice -> air"
MAT(4) = MAT(1);
MAT(4).phase  =4;

MAT(5) = MAT(2);

%
eta_min = 1e7;
eta_max = 1e30;

% ADDITIONAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAMS.YNElast=1; % elasticity on (1) or off (0)
PARAMS.YNPlas=plasYN; % plasticity on (1) or off (0)
PARAMS.tau_heal=1e12; % healing time for plasticity (s)
PARAMS.gx=0; % gravity along x
PARAMS.gy=g; % gravity along y
PARAMS.fracCFL=0.5; % distance by which a marker is allowed to move over a time step, expressed as a fraction of the smallest cell size
PARAMS.R=8.314; % gas constant
PARAMS.etamax=eta_max; % maximum viscosity
PARAMS.etamin=eta_min; % minimum viscosity
PARAMS.Tsolve=0; % yes (1) or no (0) solve for temperature
% initial temperature profile, polynomial with depth 
% T = a0 + a1*y+a2*y^2+a3*y^3+amp*sin(2*pi*X/lam)
% (make sure it matches the BCs)
PARAMS.a0=0;
PARAMS.a1=0;
PARAMS.a2=0;
PARAMS.a3=0;
PARAMS.amp=0; % amplitude of sinusoidal perturbation
PARAMS.lam=1; % wavelength of sinusoidal perturbation
PARAMS.ynTreset=1; % if ==1, reset T=T0 where im==1 (sticky layer)
PARAMS.T0=0;
% reference values for the constant diffusivity thermal solver
% (kappa = kref / (rhoref*cpref))
PARAMS.rhoref=MAT(2).rho0; 
PARAMS.kref=3;
PARAMS.cpref=1000;

% TOPOGRAPHY EVOLUTION (interface between rock and sticky air/water layer)
PARAMS.Ntopo_markers=1000; % number of markers in marker chain tracking topography
PARAMS.YNSurfaceProcesses=0; % surface processes (diffusion of topography) on or off
PARAMS.topo_kappa=1e-8; % diffusivity of topography (m^2/s)


% Solver iterations
PARAMS.Npicard_min=10; % minimum number of Picard iterations per time step
PARAMS.Npicard_max=100; % maximum number of Picard iterations per time step
PARAMS.conv_crit_ResL2=1e-9;
PARAMS.pitswitch=0; % number of Picard iterations at which the solver switches to quasi-Newton


% BOUNDARY CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pressure
PARAMS.p0cell=0; % pressure in the top-left corner of the domain (anchor point)


% flow

% boundary conditions
% entries in BC correspond to
% 1/ rollers? (1=yes, 0=no)
% 2/ type of velocity normal to boundary (0=constant)
% 3/ value of normal velocity 

BC.top=[1 0 0];
BC.bot=[BC_bot 0 0];
BC.left=[0 0 0];
BC.right=[0 0 0];

if BC.bot(1)==0
    BC.bot_profile = zeros([1 nx]);
end

if BC.top(1)==0
    BC.top_profile = zeros([1 nx]);
end


PARAMS.BalanceStickyLayer=0; % if set to 1, the code will reset the inflow 
% / outflow BCs to balance the inflow / outflow of sticky layer material,
% and rock separately, based on the position of the sticky layer / air
% interface


% thermal 

% entries in BCtherm correspond to
% 1/ type? (1=Dirichlet, 2=Neumann)
% 2/ value
BCtherm.top=[1 0];
BCtherm.bot=[1 1000];
BCtherm.left=[2 0];
BCtherm.right=[2 0];
