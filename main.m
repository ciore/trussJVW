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
statdet=checkStaticallyDeterminate(model);
if not(statdet)
  warning('The problem is statically indeterminate!')
end
figure(gcf), clf, hold on, axis equal
plotModel(model)
 
%% solve for forces
[force,react]=computeMJoints(model);
plotForce(model,force)
displayForce(react,force)

%% solve for displacements
delta=computeMVirtualWork(model,force);
plotDelta(model,delta)
displayDelta(delta)
mass=sum(model.A.*model.L.*model.rho);
deltamax=max(sqrt(sum(delta.^2,2)));



%% FUNCTIONS

%%
function plotModel(model)
  n=size(model.node,1);
  m=size(model.member,1);
  plot(model.node(:,1),model.node(:,2),'ko','MarkerFaceColor','black');
  plot(reshape(model.node(model.member,1),m,2)',reshape(model.node(model.member,2),m,2)','k')
  quiver(model.node(:,1),model.node(:,2),model.react(:,1),zeros(n,1),0.3,'g');
  quiver(model.node(:,1),model.node(:,2),zeros(n,1),model.react(:,2),0.3,'g');
  quiver(model.node(:,1),model.node(:,2),model.load(:,1),zeros(n,1),0.3,'c');
  quiver(model.node(:,1),model.node(:,2),zeros(n,1),model.load(:,2),0.3,'c');
end

%%
function plotForce(model,force)
  plot(model.node(:,1),model.node(:,2),'ko');
  h=plot(reshape(model.node(model.member(force>0,:),1),numel(find(force>0)),2)',reshape(model.node(model.member(force>0,:),2),numel(find(force>0)),2)','r');
  w=ceil(model.A(force>0,:)./max(model.A)*10);
  for i=1:numel(h)
    set(h(i),'LineWidth',w(i));
  end
  h=plot(reshape(model.node(model.member(force<0,:),1),numel(find(force<0)),2)',reshape(model.node(model.member(force<0,:),2),numel(find(force<0)),2)','b');
  w=ceil(model.A(force<0,:)./max(model.A)*10);
  for i=1:numel(h)
    set(h(i),'LineWidth',w(i));
  end
end

%%
function plotDelta(model,delta)
m=size(model.member,1);
plot(model.node(:,1)+delta(:,1),model.node(:,2)+delta(:,2),'ko');
plot(reshape(model.node(model.member,1)+delta(model.member,1),m,2)',reshape(model.node(model.member,2)+delta(model.member,2),m,2)','--k')
end

%%
function displayForce(react,force)
  fprintf(['Reaction Forces (nonzero):\n'])
  fprintf(['   node      Fx [N]\n'])
  fprintf('%7i%12.3e\n',[find(react(:,1)) nonzeros(react(:,1))]')
  fprintf(['   node      Fy [N]\n'])
  fprintf('%7i%12.3e\n',[find(react(:,2)) nonzeros(react(:,2))]')
  if exist('force')
    fprintf(['Member Forces:\n'])
    fprintf([' member       F [N]\n'])
    fprintf('%7i%12.3e\n',[(1:size(force,1)); force'])
  end
end

%%
function displayDelta(delta)
  fprintf(['Node Displacements:\n'])
  fprintf(['   node      dx [m]      dy [m]\n'])
  fprintf('%7i%12.3e%12.3e\n',[(1:size(delta,1)); delta'])
end

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
  model.L=sqrt((model.node(model.member(:,2),2)-model.node(model.member(:,1),2)).^2+(model.node(model.member(:,2),1)-model.node(model.member(:,1),1)).^2); %size mx1
  model.E=[2.1e11]*ones(size(model.member,1),1); %size mx1
  model.nu=[0.285]*ones(size(model.member,1),1); %size mx1;
  model.rho=[7800]*ones(size(model.member,1),1); %size mx1
  model.A=[0.01^2*pi/4]*ones(size(model.member,1),1); %size mx1
end

% %%
% function model=truss_def
%   model.node=[0 0; 1 0; 0.25 -0.05; 0.5 -0.05; 0.75 -0.05; 0.25 0; 0.5 0; 0.75 0]; %size nx2
%   model.member=[1 3; 3 4; 4 5; 5 2; 1 6; 6 7; 7 8; 8 2; 3 6; 3 7; 7 4; 7 5; 5 8]; %size mx2
%   model.node=[0 0; 1 0; 0.25 -0.05; 0.75 -0.05; 0.5 0]; %size nx2
%   model.member=[1 3; 3 4; 4 2; 1 5; 5 2; 3 5; 5 4]; %size mx2
%   model.react=logical([1 1; 0 1; 0 0; 0 0; 0 0; ]); %size nx2
%   model.load=[0 0; 0 0; 0 0; 0 0; 0 -1e4]; %size nx2
%   model.virt=logical([0 0; 1 0; 1 1; 1 1; 1 1]); %size nx2
%   model.L=sqrt((model.node(model.member(:,2),2)-model.node(model.member(:,1),2)).^2+(model.node(model.member(:,2),1)-model.node(model.member(:,1),1)).^2); %size mx1
%   model.E=[2.1e11]*ones(size(model.member,1),1); %size mx1
%   model.nu=[0.285]*ones(size(model.member,1),1); %size mx1;
%   model.rho=[7800]*ones(size(model.member,1),1); %size mx1
%   model.A=[0.01^2*pi/4]*ones(size(model.member,1),1); %size mx1
% end
