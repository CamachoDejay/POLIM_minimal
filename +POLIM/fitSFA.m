function [out] = fitSFA(P)
% FITSFA fits the single funnel approximation model to the input portratit.
% fit is based in a bounded fminsearch

% get modulation depths and LS from the portrait using object methods
P.getModulations;

% change 2D portrait into a its vector form
[exAngVector, emAngVector, intensity] = P.linearize;
% normalization of the intensity to its max value, probably not needed in
% the current implementation
intensity = intensity ./ max(intensity);

%%%% inputs for the fit routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   a0: initial guess for fit, [Modulation of the funnel, relative phase of
%   the funnel, geometrical ratio]. A M_f = M_ex, and a relative phase of 0
%   mean that all absorbing dipoles participate equaly in the EET emitter,
%   that is to say energy delocalization. A geometrical ratio of 1 means
%   that the side dipoles in the 3 dipole model have same amplitude as the
%   central dipole.
a0 = [P.Mex 0 1];   % Mf, thetaf, X

% Define boundary conditions:
% Lower Boundary 
%     M_f   P_f     X
LB = [0.01 -pi/2    0];
% Upper Boundary
%     M_f   P_f     X
UB = [1     pi/2    2*(1+P.Mex)/(1-P.Mex)]; 

% immutable data used whitin the minimization function. Basically data used
% to calculate RMSD - value to min.
ExpInput.Mex= P.Mex;
ExpInput.Pex = P.Pex;
ExpInput.ExAng = exAngVector;
ExpInput.EmAng = emAngVector;
ExpInput.Ftot = intensity;

% optimization settings 
opt = optimset('MaxFunEvals',10000,'MaxIter',4000,...
               'TolFun',1.0e-12,'Display','off'); % old value of tol 1e-12
% parameters not used from the fmincon function
A = []; b = []; Aeq = []; beq = []; nonlcon = [];
%%%% end of inputs for the fit routine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Actual fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[a, RMSD] = fmincon(@(a) fun2min2D(a,ExpInput),a0,...
                                            A,b,Aeq,beq,LB,UB,nonlcon,opt);
% getting output from the fit
%   modultion of the funnel
out.Mf = a(1);
%   relative phase of the funnel
out.Pf = a(2);
%   geometrical ratio used in the 3 dipole model
out.X = a(3);
%   RMSD
out.RMSD = RMSD;

% However, due to the way I constructed the fit I need to run least equares
% one last time to get the epsilon value, as in principle once you know the
% properties of the SFA there is a single answer for epsilon
model = POLIM.SFAmodel(P.Mex, P.Pex, out.Mf, out.Pf, out.X, [],...
                                                 exAngVector, emAngVector);
[ ~, out.epsilon, ~, ~ ] = SFA_fit_lsqnonneg( model, ExpInput.Ftot(:));

% generate fitted function using full range of angles 
model = POLIM.SFAmodel(P.Mex, P.Pex, out.Mf, out.Pf, out.X, out.epsilon,...
                                                 exAngVector, emAngVector);
% change from vector form to portrait (1D => 2D)
[fitPlot,exAngRad,emAngRad] = model.getPortrait;
% store data in portrait object
Pfit = POLIM.portrait(fitPlot, exAngRad, emAngRad);
% get modulation depths for the fit
Pfit.getModulations;
% store output
out.Pfit = Pfit;
        