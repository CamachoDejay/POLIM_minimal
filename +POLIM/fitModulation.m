function [Y, best_curve]=fitModulation(Angles,data2Fit,phase_res)
% FITMODULATION obtains the modulation depth, and modulation phase of the
% input curves. The input curve is assumed to have a cosine^2 dependence.
% Phase is returned in Radians!

    assert(iscolumn(Angles), 'Input Angles must be a single column vector')
    assert(size(Angles,1)==size(data2Fit,1), 'Intensity first dimentions must correspond to angles')
    
    tmp = and(all(Angles>=0), all(Angles<=pi));
    assert(tmp, 'I put angles must be in radians')
    
    m   = length(Angles);
    shapes=ones(m,2);

    if nargin == 2
        p_res = 1;

    elseif nargin == 3
        assert(and(phase_res>0,phase_res<=pi/180), 'angular resolution should be < 0.17 [rad] or 10 [degrees]')
        p_res = phase_res;
    end
    
    phases = -pi/2:p_res:pi/2;
    % as for cos^2 -pi/2 = pi/2, we remove duplicated angles.
    phases(phases == pi/2) = [];

    % unit important variables for the loop
    min_error  = inf(1,size(data2Fit,2));
    best_curve = NaN(size(data2Fit));
    best_m = inf(1,size(data2Fit,2));
    best_phase = NaN(1,size(data2Fit,2));
    %%Fitting
     for i = 1:length(phases)
         % modify model for linear fit, 1st col is always flat background,
         % 2nd col is the cos^2 at certain angular phase in rad
        phase       = phases(i);
        shapes(:,2) = cos( 2 * (Angles - phase));
        % fit
        x       = lscov(shapes,data2Fit);
        % extra calculations to ensure that it is a positive solution, it
        % is left this way as I can make it parallel for many input curves
        eps     = abs(x);
        modulations = eps(2,:) ./ eps(1,:);
        
        y_fit   = shapes*eps;
        y_error = data2Fit-y_fit;
        y_error = sum(y_error.^2,1);
        
        % see if this fit is better than previous ones
        idx = y_error < min_error;
        % update ouput accordingly
        min_error(idx) = y_error(idx);
        best_curve(:,idx) = y_fit (:,idx);
        best_m(idx) = modulations(idx);
        best_phase(idx) = phases(i);

    end
    %%Fitting output    
    Y = cat(1, best_m, best_phase, min_error);
    
end