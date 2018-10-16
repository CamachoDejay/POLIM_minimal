clear
close all
clc

% let us first load data coming from homo-FRET between EGFP molecules in a
% cubic latice as calculated by:
% https://github.com/CamachoDejay/FRET-calculations and presented in the
% article: "2D polarization imaging as a low-cost fluorescence method to
% detect ?-synuclein aggregation ex vivo in models of Parkinson?s disease"
% for simplicity I have uploaded 4 portraits to the data folder

% load portraits, excitation and emission angles
load([cd filesep 'data' filesep 'portraits.mat'])
load([cd filesep 'data' filesep 'emAngRad.mat'])
load([cd filesep 'data' filesep 'exAngRad.mat'])

% lets now iterate over the portraits, the interchromophoric
% distance decreases with increasing index, therefore epsilon should also
% decrease as it will measure the homo-FRET efficiency.
for i = 1:4
    % get a particular portrait data, 1-dim: excitation angles index,
    % 2-dim: emission angles index
    Intensity = portraits(:,:,i);
    % build a portrait object, that contains the data and can do some basic
    % calculations, such as modulation depths
    P = POLIM.portrait(Intensity, exAngRad, emAngRad);
    % get modulation depths
    P.getModulations;
    % show the portrait
    P.showPortrait;
    % modify figure's title
    a = gca;
    titleStr = sprintf('Experimental portrait, idx: %.0f', i);
    mStr = sprintf('M_e_x: %.2f, M_e_m: %.2f, LS: %.2f',P.Mex, P.Mem, P.LS);
    epsStr = sprintf('Epsilon: To be determined');
    a.Title.String = {titleStr, mStr, epsStr};

    % now we run the SFA fitting on the portrait
    [SFAoutput] = POLIM.fitSFA(P);
    % the fitted portrait is stored so we can also show it easily
    SFAoutput.Pfit.showPortrait;
    % lets also modify the figure' title so it is easy to indentify
    a = gca;
    titleStr = sprintf('Fitted portrait, idx: %.0f', i);
    mStr = sprintf('M_e_x: %.2f, M_e_m: %.2f, LS: %.2f',...
        SFAoutput.Pfit.Mex, SFAoutput.Pfit.Mem, SFAoutput.Pfit.LS);
    epsStr = sprintf('Epsilon: %.2f',SFAoutput.epsilon);
    a.Title.String = {titleStr, mStr, epsStr};

end

%% as a second example let us create a portrait under the SFA and try to 
% recover the input values of modulation and energy transfer.

% generating external SFA model
%   number of excitaiton angles
nExAng = 180;
%   number of emission angles
nEmAng = 180;
%   modulation in excitation
Mex = .2;
%   modulation phase in excitation
Pex = 0*pi/180;
%   modulation of the funnel
Mf  = .4;
%   phase of the funnel relative to the excitation
Pf  = 0*pi/180;
%   geometrical ratio
X = (1+Mex)/(1-Mex);
%   epsilon - EET metric
epsilon = .5;
%   creating empty portrait object
Intensity = zeros(nExAng,nEmAng);
exAngRad  = linspace(0,pi,nExAng+1);
exAngRad(end) = [];
emAngRad  = linspace(0,pi,nEmAng+1);
emAngRad(end) = [];
P = POLIM.portrait(Intensity, exAngRad, emAngRad);
%   get vector version of the portrait
[exAngVector, emAngVector, ~] = P.linearize;
%   using the right angles, get vector version of the SFA
model = POLIM.SFAmodel(Mex, Pex, Mf, Pf, X, epsilon, exAngVector, emAngVector);
[fitPlot,exAngRad,emAngRad] = model.getPortrait;
% store data back in portrait object
P = POLIM.portrait(fitPlot, exAngRad, emAngRad);

% now we fit the external portrait with the SFA fit function
%   get modulation depths
P.getModulations;
%   show portrait
P.showPortrait;
%   modify title of the figure
a = gca;
mStr = sprintf('M_e_x: %.2f, M_e_m: %.2f, LS: %.2f',P.Mex, P.Mem, P.LS);
epsStr = sprintf('Epsilon: To be determined');
a.Title.String = {'Simulated portrait', mStr, epsStr};
%   do SFA fit
[SFAoutput] = POLIM.fitSFA(P);
%   show output
SFAoutput.Pfit.showPortrait;
%   modify title
a = gca;
mStr = sprintf('M_e_x: %.2f, M_e_m: %.2f, LS: %.2f',...
    SFAoutput.Pfit.Mex, SFAoutput.Pfit.Mem, SFAoutput.Pfit.LS);
epsStr = sprintf('Epsilon: %.2f',SFAoutput.epsilon);
a.Title.String = {'Fited portrait', mStr, epsStr};
