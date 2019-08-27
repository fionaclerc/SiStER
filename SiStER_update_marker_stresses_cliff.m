%=============================================================================
% modified from "SiStER_update_marker_stresses" to track topography
% of rectangular cliff (F. Clerc 2018)
%
% Updates stresses on markers for CURRENT solution.  Stress rotation
% occurs after solutions are output
% G.Ito 8/16
%==============================================================================

% Compute STRESS Changes on nodes, interpolate to markers, and apply to marker stresses
dsxx=(2*etan.*EXX-sxxOLD).*Zn;
dsxy=(2*etas.*EXY-sxyOLD).*Zs;

[dsxxm]=SiStER_interp_normal_nodes_to_markers(dsxx,xc,yc,xm,ym,icn,jcn);
[dsxym]=SiStER_interp_shear_nodes_to_markers(dsxy,x,y,xm,ym,icn,jcn);
sxxm=sxxm+dsxxm;
sxym=sxym+dsxym;

% topography markers
[dsxxm1]=SiStER_interp_normal_nodes_to_markers(dsxx,xc,yc,topo_x1,topo_y1,icn1,jcn1);
[dsxym1]=SiStER_interp_shear_nodes_to_markers(dsxy,x,y,topo_x1,topo_y1,icn1,jcn1);
sxxm1=sxxm1+dsxxm1;
sxym1=sxym1+dsxym1;

[dsxxm2]=SiStER_interp_normal_nodes_to_markers(dsxx,xc,yc,topo_x2,topo_y2,icn2,jcn2);
[dsxym2]=SiStER_interp_shear_nodes_to_markers(dsxy,x,y,topo_x2,topo_y2,icn2,jcn2);
sxxm2=sxxm2+dsxxm2;
sxym2=sxym2+dsxym2;

