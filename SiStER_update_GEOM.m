function [im,GEOM] = SiStER_update_GEOM(vel,dt,GEOM,xm,ym,Nphase,direction,h,H)
% change geometry at rate "vel" in m/s over timestep dt (F. Clerc 2018)
% make sure timesteps are large enough relative to grid size
% direction = 1 if going up, direction = -1 if going down.
% direction = 0 if both.

if direction == -1
    GEOM(5).top = min(GEOM(5).top + (vel.*dt), GEOM(5).bot);
elseif direction == 1
    GEOM(5).bot = max(GEOM(5).bot - (vel.*dt), GEOM(5).top);
elseif direction == 0
    GEOM(5).bot = max(GEOM(5).bot - (vel.*(1-(h./H)).*dt), GEOM(4).bot);
    GEOM(5).top = min(GEOM(5).top + (vel.*((h./H)).*dt), GEOM(4).bot);
end
[im] = SiStER_initialize_marker_phases(Nphase,GEOM,xm,ym);