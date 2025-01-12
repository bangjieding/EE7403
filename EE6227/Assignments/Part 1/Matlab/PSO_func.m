function [gbest, gbestval, fitcount, gbestval_con]= PSO_func(fhd, Dimension, Particle_Number, Max_Gen, VRmin, VRmax, Inertia_Weight, C_set, varargin)
    me = Max_Gen;
    ps = Particle_Number;
    D = Dimension;
    iwt = Inertia_Weight;
    cc  = C_set;
    % cc=[2 2];   %acceleration constants
    % iwt=0.9-(1:me).*(0.5./me); %inertia fector table
    % iwt=0.5.*ones(1,me);
    if length(VRmin) == 1
        VRmin = repmat(VRmin,1,D);
        VRmax = repmat(VRmax,1,D);
    end
    mv = 0.2 * (VRmax - VRmin);
    VRmin = repmat(VRmin, ps, 1);
    VRmax = repmat(VRmax, ps, 1);
    Vmin = repmat(-mv, ps, 1);
    Vmax = -Vmin;
    pos = VRmin + (VRmax - VRmin) .* rand(ps, D);

    e = feval(fhd, pos', varargin{:});

    fitcount = ps;
    vel = Vmin+2.*Vmax.*rand(ps, D);%initialize the velocity of the particles
    pbest = pos;
    pbestval = e; %initialize the pbest and the pbest's fitness value
    [gbestval, gbestid] = min(pbestval);
    gbest = pbest(gbestid, :);%initialize the gbest and the gbest's fitness value
    gbestrep = repmat(gbest, ps, 1);
    gbestval_con(1) = gbestval;
    for i = 2:me
            aa = cc(1) .* rand(ps, D) .* (pbest - pos) + cc(2) .* rand(ps, D) .* (gbestrep - pos);
            vel = iwt(i) .* vel + aa;
            vel = (vel > Vmax) .* Vmax + (vel <= Vmax) .* vel;
            vel = (vel < Vmin) .* Vmin + (vel >= Vmin) .* vel;
            pos = pos + vel;
            pos = ((pos >= VRmin) & (pos <= VRmax)) .* pos...
                + (pos < VRmin) .* (VRmin + 0.25 .* (VRmax - VRmin) .* rand(ps, D)) + (pos > VRmax) .* (VRmax - 0.25 .* (VRmax - VRmin) .* rand(ps, D));
            e = feval(fhd, pos', varargin{:});
            fitcount = fitcount + ps;
            tmp = (pbestval < e);
            temp = repmat(tmp', 1, D);
            pbest = temp .* pbest + (1 - temp) .* pos;
            pbestval = tmp .* pbestval + (1 - tmp) .* e;%update the pbest
            [gbestval, tmp] = min(pbestval);
            gbest = pbest(tmp, :);
            gbestrep = repmat(gbest, ps, 1);%update the gbest
            gbestval_con(i) = gbestval;
        end
    end


