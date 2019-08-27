% modified from "SiStER_update_topography_markers"
% updated to keep track of cliff topography (F. Clerc 2018)

% advect the marker chain that keeps track of topography 
% in the current flow field
[topo_x1,topo_y1] = SiStER_advect_markers(x,y,topo_x1,topo_y1,dx,dy,dt_m,vx,vy);
[topo_x2,topo_y2] = SiStER_advect_markers(x,y,topo_x2,topo_y2,dx,dy,dt_m,vx,vy);

% locate the interface between sticky layer and left / right edge
% 
if isempty(find(topo_x1<0,1))==1
    topoL=topo_y1(1);
else
    topoL=interp1(topo_x1,topo_y1,0);
end

if isempty(find(topo_y2>ysize,1))==1
    topoR=topo_x2(end);
else
    topoR=interp1(topo_y2,topo_x2,ysize);% change this
end

% eliminate topography markers that left domain, keep the first one out on both sides
Iin1=find(topo_x1>0);
Iin2=find(topo_y2<ysize);

topo_x1=topo_x1(Iin1);
topo_y1=topo_y1(Iin1);
topo_x2=topo_x2(Iin2);
topo_y2=topo_y2(Iin2);

topo_x1=[0 topo_x1];
topo_y1=[topoL topo_y1];
topo_x2 = [topo_x2 topoR];
topo_y2 = [topo_y2 ysize];

if PARAMS.YNSurfaceProcesses==1
    % ERODE TOPOGRAPHY
    [topo_y]=SiStER_topography_diffusion_solver(topo_x,topo_y,dt_m,PARAMS.topo_kappa);
    % RESET ROCK AND AIR (assumes topography is only interface between phase 1 and 2)
    topomarkers=interp1(topo_x,topo_y,xm);
    im(im==1 & ym>=topomarkers)=2;
    im(im>=2 & ym<topomarkers)=1;
end

% if there has been too much stretching, regrid the surface topography
if max(diff(topo_x1))>5*topo_marker_spacing1 || issorted(topo_x1)==0
    % surface regridding happens if somewhere 2 topo markers have been
    % stretched apart by more than 5 times the inital mean marker spacing
    % or if topo_x is no longer sorted due to compression.
    topo_xREGRID = linspace(0,GEOM(3).left,Ntopo/2);
    topo_yREGRID=interp1(topo_x1,topo_y1,topo_xREGRID(2:end));
    topo_yREGRID=[topoL topo_yREGRID];
    topo_x1=topo_xREGRID;
    topo_y1=topo_yREGRID;
    disp('**REGRIDDING TOPOGRAPHY MARKERS**')
end

% if there has been too much stretching, regrid the surface topography
if max(diff(topo_y2))>5*topo_marker_spacing2 || issorted(topo_y2)==0
    % surface regridding happens if somewhere 2 topo markers have been
    % stretched apart by more than 5 times the inital mean marker spacing
    % or if topo_x is no longer sorted due to compression.
    topo_yREGRID = linspace(GEOM(2).top,GEOM(2).bot,Ntopo/2);
    topo_xREGRID=interp1(topo_y2,topo_x2,topo_yREGRID(1:end-1));
    topo_xREGRID=[topo_xREGRID topoR];
    topo_x2=topo_xREGRID;
    topo_y2=topo_yREGRID;
    disp('**REGRIDDING TOPOGRAPHY MARKERS**')
end


[qd1,icn1,jcn1] = SiStER_locate_markers_in_grid(topo_x1,topo_y1,x,y,dx,dy);
[qd2,icn2,jcn2] = SiStER_locate_markers_in_grid(topo_x2,topo_y2,x,y,dx,dy);
topo_x = [topo_x1 topo_x2];topo_y = [topo_y1 topo_y2];


