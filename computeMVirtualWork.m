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
% computeMVirtualWork.m
% Compute dispacement of nodes using method of virtual work
% Must precompute forces in members first

function delta=computeMVirtualWork(model,force)
  n=size(model.node,1);
  delta=zeros(n,2);
  length=sqrt((model.node(model.member(:,2),2)-model.node(model.member(:,1),2)).^2+(model.node(model.member(:,2),1)-model.node(model.member(:,1),1)).^2);
  for i=1:2*n
    if model.virt(i)
      model.load(:)=0;
      model.load(i)=1;
      virtual=computeMJoints(model);
      delta(i)=sum(virtual.*force.*length./model.A./model.E);
    end
  end
end
