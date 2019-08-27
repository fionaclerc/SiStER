% SiStER_MAIN_cliff.m
%
% Reproduces results of :
% F. Clerc, B.M. Minchew, M.D. Behn (2019). Marine Ice Cliff Instability
%   Mitigated by Slow Removal of Ice Shelves.
% fclerc <at> mit.edu
% Aug. 2019
%
% Inputs :
%   1. "fname" - name of output files.
%   2. "deltat" - timescale of ice-shelf removal
%   3. "h" - subaerial cliff height
%   4. "n_flow" - stress exponent in flow law (use 1 or 3)
%
%
%
% modified from SiStER_MAIN.m :
%
%   Simple Stokes solver with Exotic Rheologies
%
%   Main routine doing initialization, time loop and outputs
%
%   J.-A. Olive, B.Z. Klein, E. Mittelstaedt, M. Behn, G. Ito, S. Howell
%   jaolive <at> ldeo.columbia.edu
%   March 2011 - April 2017
%
%


function SiStER_MAIN_cliff(fname,deltat,h,n_flow)



BC_bot = 1;% bottom boundary condition; 0: no slip, 1: free slip

Dt = 1e3;% non-dimensionalized time scale
nt_trans = 50;% number of timesteps over transition
dt_m2 = Dt./nt_trans;% min. step size after initiation of transition [s]
dt_m0 = 100;% min. step size before transition [s]

Tfinal = 1.2e3;% end time [s]
trans = 2;
G0 = 2e9;% dimensionalized shear modulus [Pa]
Gice = G0.*deltat./Dt;% non-dimensionalized shear modulus

direction = 0;
dxx = 20;
dyy = h./nt_trans;% must make vertical cells small enough such that
% ice-shelf is thinned linearly (not step-wise).

% plasticity is turned off.
tau_yield = 1e12;plasYN = 0;


InpFil = 'SiStER_Input_File_coulomb_cliff_GEOM.m';

g = 9.8;% gravitational acceleartion [m/s^2]
rhoi = 900;rhow = 1000;% density of ice/water [kg/m^3]
H = h./(1 - rhoi./rhow); % Total ice thickness.


A0 = 1.2e-25;% Glen's flow-law parameter.

if n_flow==3
    A = A0;
elseif n_flow==1
    % average deviatoric stress
    tzz = @(h,H) ((g./4).*(H.*(rhoi-rhow) + h.*rhow.*(2 - h./H)));
    
    % effective viscosity
    eta_eff2 = @(te,A) te.^-2./(2.*A);
    
    % modify pre-factor to include effective stress term
    A = 1./(2.*eta_eff2(tzz(h,H),A0));
end

RAS = 0.01;% Tensile/Coulombic threshold R.


t_trans = [dt_m0 dt_m0+Dt];

% INITIALIZATION

dt_m= dt_m0;

% Input File: loads parameter values, model geometry, boundary conditions
run(InpFil)

% construct grid and initialize marker / node arrays
SiStER_Initialize_cliff_topo;

% BEGIN TIME LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time = -dt_m0;

vel = H./diff(t_trans);% rate of thinning of buttressing ice-shelf.
ver = 6;

if ~isfield('BC','midx')
    BC.midx = [0 0 0];
end

%%
Nt = ceil(Tfinal./dt_m2);

for t=[1 1:Nt*1e3] % time loop
    if time<Tfinal
    
    
    disp(['STARTING ITERATION: ' num2str(t) ' out of ' num2str(Nt) sprintf(' time: %.2e of %.2e',time,Tfinal)])
    
    % update time
    time=time+dt_m;
    
     % Here we prepare nodal arrays to feed the Stokes solver 
    SiStER_material_props_on_nodes

    %%% SOLVE STOKES WITH NON-LINEAR RHEOLOGY HERE 
    %SiStER_flow_solve_sigma
    SiStER_flow_solve
    
    
    % GET STRAIN RATE FROM CURRENT SOLUTION
    epsIIm=SiStER_interp_shear_nodes_to_markers(epsII_s,x,y,xm,ym,icn,jcn);
    epsIIm1=SiStER_interp_shear_nodes_to_markers(epsII_s,x,y,topo_x1,topo_y1,icn1,jcn1);
    epsIIm2=SiStER_interp_shear_nodes_to_markers(epsII_s,x,y,topo_x2,topo_y2,icn2,jcn2);
    
    % USE STRAIN RATE TO UPDATE STRESSES ON MARKERS
    SiStER_update_marker_stresses_cliff;
    
    % BUILD UP PLASTIC STRAIN IN YIELDING AREAS IF PLASTICITY IS ACTIVATED
    if (PARAMS.YNPlas==1) 
        SiStER_update_ep;
    end
  
    % OUTPUT VARIABLES OF INTEREST (prior to rotation & advection)
    if (mod(t,dt_out)==0 && dt_out>0) || t==0 || t==Nt % SAVING SELECTED OUTPUT
        disp('SAVING SELECTED VARIABLES TO OUTPUT FILE') 
        filename = sprintf('%03d%s',t,fname);
        [etam]=SiStER_interp_shear_nodes_to_markers(etas,x,y,xm,ym,icn,jcn); % to vsisualize viscosity on markers
        save(filename,'X','Y','vx','vy','p','time','xm','ym','etam',...
            'rhom','BC','im','epsIIm','epsIIm1','epsIIm2','epsII_s',...
            'sxxm','sxym','sxxm1','sxym1','sxxm2','sxym2','rho',...
            'sxxOLD','sII','topo_x1','topo_y1','topo_x2','topo_y2',...
            'topo_x','topo_y','ResL2');

    end
    
    
    % SET ADVECTION TIME STEP BASED ON CURRENT FLOW SOLUTION
    [dt_m]=SiStER_set_timestep(dx,dy,vx,vy,PARAMS);
    
    
    % USE PRE-SET MIN. TIME STEPS,
    if time<t_trans(1)
        dt_m = dt_m0;
    elseif dt_m2>0
        if dt_m < dt_m2
            disp(sprintf('dt_m is small : %.1e',dt_m));
        else
            dt_m =dt_m2;
        end
    end
    
    % ROTATE ELASTIC STRESSES IN CURRENT FLOW FIELD
    if (PARAMS.YNElast==1) 
        SiStER_rotate_stresses;
    end
    
    % EVOLVE TEMPERATURE FIELD THROUGH DIFFUSION
    if PARAMS.Tsolve==1
        SiStER_thermal_update;
    end

    % MARKER ADVECTION, REMOVAL, AND ADDITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SiStER_move_remove_and_reseed_markers;
    % advect markers in current flow field
    % remove markers if necessary
    % add markers if necessary
    SiStER_update_topography_markers_cliff;
    % here we do the same for the marker chain that keeps track of topography
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    % THIN THE BUTTRESSING ICE-SHELF LINEARLY.
    if time>= t_trans(1)
        im = SiStER_update_GEOM(vel,time - t_trans(1),GEOM,xm,ym,Nphase,direction,h,H);
    end

    
    disp('---------------')
    disp(['END OF ITERATION: ' num2str(t) ' out of ' num2str(Nt) ' - SIMULATION TIME: ' num2str(time/365.25/24/3600/1000) ' kyrs.'])
    disp('--------------------------------')
    disp('--------------------------------')
    

    end
        
end

disp('FIN');

    