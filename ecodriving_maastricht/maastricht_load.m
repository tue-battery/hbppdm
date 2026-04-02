function [mstr_data] = maastricht_load(Ts)

    load Verification_route.mat;
    
    distance_route_m = Distance_fit*1000;
    length_route = Distance_fit(end)*1000;
    Distance_interp = 0:Ts:length_route;
    route_len = length(Distance_interp);

    Height_interp = interp1(distance_route_m,Height_fit_smooth,Distance_interp,"linear","extrap");
    vmax_interp   = interp1(distance_route_m,v_max,Distance_interp,"linear","extrap");

    mstr_data.distance = Distance_interp;
    mstr_data.height   = Height_interp;
    mstr_data.vmax     = vmax_interp;
    mstr_data.N        = route_len-1;
    mstr_data.Ts       = Ts;

end