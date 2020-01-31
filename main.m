% This file is part of trussJVW, a code to compute the forces and
% deflections in a truss using the method of joints and the method
% virtual work
% 
% Copyright (C) 2020 Ciarán O'Reilly <ciaran@kth.se>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% main.m 

clear

%% define truss model and visualise
model=truss_def; %call function defining truss problem
statdet=trussJVWFuncs.checkStaticallyDeterminate(model);
figure(gcf), clf, hold on, axis equal
trussJVWFuncs.plotModel(model)
 
%% solve for forces
[force,react]=trussJVWFuncs.computeMJoints(model);
trussJVWFuncs.plotForce(model,force)
trussJVWFuncs.displayForce(react,force)

%% solve for displacements
delta=trussJVWFuncs.computeMVirtualWork(model,force);
trussJVWFuncs.plotDelta(model,delta)
trussJVWFuncs.displayDelta(delta)
mass=sum(model.A.*model.L.*model.rho);
deltamax=max(sqrt(sum(delta.^2,2)));



%% FUNCTIONS

%%
function model=tie_def
  model.node=[0 0; 1 0]; %size nx2
  model.member=[1 2]; %size mx2
  model.react=logical([1 1; 0 1]); %size nx2
  model.load=[0 0; 1e4 0]; %size nx2
  model.virt=logical([0 0; 1 0]); %size nx2
  model.L=sqrt((model.node(model.member(:,2),2)-model.node(model.member(:,1),2)).^2+(model.node(model.member(:,2),1)-model.node(model.member(:,1),1)).^2);  %size mx1
  model.E=[2.1e11]*ones(size(model.member,1),1); %size mx1
  model.nui=[0.285]*ones(size(model.member,1),1); %size mx1;
  model.rho=[7800]*ones(size(model.member,1),1); %size mx1
  model.A=[0.01^2*pi/4];  %size mx2
end

%%
function model=truss_def
  model.node=[0 0; 1 0; 0.25 -0.05; 0.5 -0.05; 0.75 -0.05; 0.25 0; 0.5 0; 0.75 0]; %size nx2
  model.member=[1 3; 3 4; 4 5; 5 2; 1 6; 6 7; 7 8; 8 2; 3 6; 3 7; 7 4; 7 5; 5 8]; %size mx2
  model.react=logical([1 1; 0 1; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0]); %size nx2
  model.load=[0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 -1e4; 0 0]; %size nx2
  model.virt=logical([0 0; 1 0; 1 1; 1 1; 1 1; 1 1; 1 1; 1 1]); %size nx2
  model.L=sqrt((model.node(model.member(:,2),2)-model.node(model.member(:,1),2)).^2+(model.node(model.member(:,2),1)-model.node(model.member(:,1),1)).^2)'; %size mx1
  model.E=[2.1e11]*ones(size(model.member,1),1); %size mx1
  model.nu=[0.285]*ones(size(model.member,1),1); %size mx1;
  model.rho=[7800]*ones(size(model.member,1),1); %size mx1
  model.A=[0.01^2*pi/4]*ones(size(model.member,1),1); %size mx1
end
