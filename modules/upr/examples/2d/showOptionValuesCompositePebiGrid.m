%% Example
% In this example we show all optional options for compositePebiGrid.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%% wellLines
% This sets the wells in our reservoir. The wells are stored as a cell
% array, each element corresponding to one well
w = {[0.2,0.8;0.5,0.6;0.8,0.8],...
     [0.5,0.2]};
gS = [0.1,0.1];
pdims=[1,1];
G = compositePebiGrid(gS, pdims,'wellLines',w);

clf
plotGrid(G);
plotGrid(G,G.cells.tag,'facecolor','b')
axis equal tight off
plotLinePath(w,'wo-','linewidth',2,'MarkerFaceColor','w');


%% wellGridFactor
% The wellGridFactor sets the relative size of the well cells. If
% wellGridFactor=0.5 the well cells will have about half the size of the
% reservoir sites:
w = {[0.2,0.3;0.8,0.7]};
gS = [0.1,0.1];
pdims=[1,1];
G1 = compositePebiGrid(gS, pdims,'wellLines',w,'wellGridFactor',1);
G2 = compositePebiGrid(gS, pdims,'wellLines',w,'wellGridFactor',1/2);

figure('Position',[480 340 980 420])
subplot(1,2,1)
plotGrid(G1), title('wellGridFactor=1'), axis equal tight off
plotGrid(G1,G1.cells.tag,'facecolor','b')
hold on, plotLinePath(w,'wo-','linewidth',2,'MarkerFaceColor','w');
subplot(1,2,2)
plotGrid(G2), title('wellGridFactor=1/2'), axis equal tight off
plotGrid(G2,G2.cells.tag,'facecolor','b')
hold on, plotLinePath(w,'wo-','linewidth',2,'MarkerFaceColor','w');


%% mlqtMaxLevel
% mlqtMaxLevel sets how many steps of local grid refinement we create
% around each well.

w = {[0.2,0.3;0.5,0.5;0.8,0.5]};
gS = [0.1,0.1];
pdims=[1,1];
G1 = compositePebiGrid(gS, pdims,'wellLines',w,'wellGridFactor',1/2,'mlqtMaxLevel',1);
G2 = compositePebiGrid(gS, pdims,'wellLines',w,'wellGridFactor',1/8,'mlqtMaxLevel',3);

figure('Position',[480 340 980 420])
subplot(1,2,1)
plotGrid(G1); title('mlqtMaxLevel=1'), axis equal tight off
plotGrid(G1,G1.cells.tag,'facecolor','b')
hold on, plotLinePath(w,'wo-','linewidth',2,'MarkerFaceColor','w');
subplot(1,2,2)
plotGrid(G2); title('mlqtMaxLevel=3'), axis equal tight off
plotGrid(G2,G2.cells.tag,'facecolor','b')
hold on, plotLinePath(w,'wo-','linewidth',2,'MarkerFaceColor','w');


%% mlqtLevelSteps
% We can define which distance from a well triggers which level of
% mlqt-refinement.
w = {[0.2,0.3;0.5,0.5;0.8,0.5]};
gS = [0.1,0.1];
pdims=[1,1];
lev1 = [0.05,0.1];
lev2 = [0.1,0.2];
G1 = compositePebiGrid(gS, pdims,'wellLines',w,'wellGridFactor',1/4, ...
                      'mlqtMaxLevel',2,'mlqtLevelSteps',lev1);
G2 = compositePebiGrid(gS, pdims,'wellLines',w,'wellGridFactor',1/4, ...
                      'mlqtMaxLevel',2,'mlqtLevelSteps',lev2);

figure('Position',[480 340 980 420])
subplot(1,2,1)
plotGrid(G1); title('mlqtLevelSteps=[0.05,0.1]'), axis equal tight off
plotGrid(G1,G1.cells.tag,'facecolor','b')
hold on, plotLinePath(w,'wo-','linewidth',2,'MarkerFaceColor','w');
subplot(1,2,2)
plotGrid(G2); title('mlqtLevelSteps=[0.1,0.2]'), axis equal tight off
plotGrid(G2,G2.cells.tag,'facecolor','b')
hold on, plotLinePath(w,'wo-','linewidth',2,'MarkerFaceColor','w');


%% wellRho
% A function that sets the relative size of the fault sites in the domain.
w = {[0.2,0.3;0.5,0.5;0.8,0.7]};
gS = [0.1,0.1];
pdims=[1,1];
wellRho = @(p) 1 - 0.9*p(:,1);
G1 = compositePebiGrid(gS, pdims,'wellLines',w);
G2 = compositePebiGrid(gS, pdims,'wellLines',w,'wellRho',wellRho);

figure('Position',[480 340 980 420])
subplot(1,2,1)
plotGrid(G1); plotGrid(G1,G1.cells.tag,'facecolor','b');
plotGrid(G1,G1.cells.tag,'facecolor','b')
hold on, plotLinePath(w,'wo-','linewidth',2,'MarkerFaceColor','w');
title('wellRho=1'); axis equal tight off
subplot(1,2,2)
plotGrid(G2); plotGrid(G2,G2.cells.tag,'facecolor','b');
plotGrid(G2,G2.cells.tag,'facecolor','b')
hold on, plotLinePath(w,'wo-','linewidth',2,'MarkerFaceColor','w');
title('wellRho=@(p) 1 - 0.9*p(:,1)')
axis equal tight off


%% protLayer
% Adds a protection layer around wells. The protection sites are placed
% normal along the well path
x = linspace(0.2,0.8);
y = 0.5+0.1*sin(pi*x);
w = {[x',y']};
gS = [0.1,0.1];
pdims=[1,1];
G1 = compositePebiGrid(gS, pdims, 'wellLines', w, ...
                       'interpolWP', true, 'protLayer', false);
G2 = compositePebiGrid(gS, pdims, 'wellLines', w, ...
                       'interpolWP', true, 'protLayer', true);

figure('Position',[480 340 980 420])
subplot(1,2,1)
plotGrid(G1); plotGrid(G1,G1.cells.tag,'facecolor','b')
title('Protection Layer off'), axis equal tight off
subplot(1,2,2)
plotGrid(G2); plotGrid(G2,G2.cells.tag,'facecolor','b')
title('Protection Layer on'), axis equal tight off


%% protD
% This parameter sets the distance from the protection sites to the well
% path. This is cellarray of functions, one for each well path.
w = {[0.2,0.8;0.8,0.8],...
     [0.5,0.2;0.1,0.5]};
gS = [0.08,0.08];
pdims=[1,1];
protD = {@(p)0.01+0.1*p(:,1), @(p) 0.01*ones(size(p,1),1)};
G1 = compositePebiGrid(gS, pdims,'wellLines',w,'protLayer',true);
G2 = compositePebiGrid(gS, pdims,'wellLines',w,'protLayer',true,'protD',protD);

figure('Position',[480 340 980 420])
subplot(1,2,1)
plotGrid(G1); plotGrid(G1,G1.cells.tag,'facecolor','b')
title('Default protection distance'), axis equal tight off
subplot(1,2,2)
plotGrid(G2); plotGrid(G2,G2.cells.tag,'facecolor','b')
title('User-prescribed protection distance'), axis equal tight off


%% faultLines
% This sets the faults in our reservoir. The faults are stored as a cell
% array, each element corresponding to one fault
f = {[0.2,0.8;0.5,0.65;0.8,0.8],...
     [0.5,0.2;0.1,0.5]};
gS = [0.1,0.1];
pdims=[1,1];
figure()
G1 = compositePebiGrid(gS, pdims,'faultLines',f);
G2 = compositePebiGrid(gS, pdims,'faultLines',f,'interpolFL',[true; true]);

figure('Position',[480 340 980 420])
subplot(1,2,1)
plotGrid(G1); plotFaces(G1,G1.faces.tag,'edgeColor','r','LineWidth',4)
plotLinePath(f,'--ow','linewidth',2,'MarkerFaceColor','w');
title('Fault lines: no interpolation'), axis equal tight off
subplot(1,2,2)
plotGrid(G2); plotFaces(G2,G2.faces.tag,'edgeColor','r','LineWidth',4)
plotLinePath(f,'--ow','linewidth',2,'MarkerFaceColor','w');
title('Fault lines: interpolated'), axis equal tight off


%% faultGridFactor
% The faultGridFactor sets the relative distance between the fault cells
% along the fault paths. If fautlGridFactor=0.5 the fault cells will be
% have spaceing about half the size of the reservoir cells:
f = {[0.2,0.3;0.5,0.5;0.8,0.5]};
gS = [0.1,0.1];
pdims=[1,1];
G1 = compositePebiGrid(gS, pdims,'faultLines',f,'faultGridFactor',1);
G2 = compositePebiGrid(gS, pdims,'faultLines',f,'faultGridFactor',1/2);

figure('Position',[480 340 980 420])
subplot(1,2,1)
plotGrid(G1); title('faultGridFactor=1'); axis equal tight off
plotLinePath(f,'--or','linewidth',2,'MarkerFaceColor','w');
subplot(1,2,2)
plotGrid(G2); title('faultGridFactor=1/2'); axis equal tight off
plotLinePath(f,'--or','linewidth',2,'MarkerFaceColor','w');


%% circleFactor
% The fault sites are generated by setting placing a set of circles along
% the fault path with a distance gS*faultGridFactor. The radius of the
% circles are gS*faultGridFactor*circleFactor. The fault sites are place
% at the intersection of these circles.
f = {[0.2,0.3;0.5,0.5;0.8,0.5]};
gS = [0.1,0.1];
pdims=[1,1];
G1 = compositePebiGrid(gS, pdims,'faultLines',f,'circleFactor',0.55);
G2 = compositePebiGrid(gS, pdims,'faultLines',f,'circleFactor',0.9);

figure('Position',[480 340 980 420])
subplot(1,2,1)
plotGrid(G1); title('circleFactor=0.55'), axis equal tight off
plotLinePath(f,'--or','linewidth',2,'MarkerFaceColor','w');
subplot(1,2,2)
plotGrid(G2); title('circleFactor=0.9'), axis equal tight off
plotLinePath(f,'--or','linewidth',2,'MarkerFaceColor','w');

%% PolyBdr
% Create a non-square reservoir domain
bdr   = [0,0;-0.5,1;0.5,1.5;1,1];
gs    = [0.1,0.1];
pdims = [1,1]; % can be set to anything
G = compositePebiGrid(gs,pdims,'polyBdr',bdr);

figure()
plotGrid(G); title('Polygonal boundary'); axis equal tight off
