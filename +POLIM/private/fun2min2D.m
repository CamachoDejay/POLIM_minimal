function RMSD = fun2min2D(a, ExpInput)
%FUN2MIN2D function used during minimization to find best fit of data. It
%creates a SFA model, find best epsilon using least squares, and returns
%the RMSD

% get values to fit, that is to say values that change with each iteration
MF    = a(1);
thetaf = a(2);
X     = a(3);

% get fixed values - data to fit
Mex = ExpInput.Mex;
Pex = ExpInput.Pex;
ExAng = ExpInput.ExAng;
EmAng = ExpInput.EmAng;

% now the SFA model to use depends on both the fit values and the fixed
% values
model = POLIM.SFAmodel(Mex, Pex, MF, thetaf, X, [], ExAng, EmAng);

% now that I have my fully constrained model, find goodnes of fit
[ ~, ~, RMSD, ~ ] = SFA_fit_lsqnonneg( model, ExpInput.Ftot(:));