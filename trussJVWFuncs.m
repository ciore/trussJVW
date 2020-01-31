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

classdef trussJVWFuncs
  methods(Static)
    
    %%
    function statdet=checkStaticallyDeterminate(model)
    % Checks if Nmembers+Nreactions==2*Nnodes
      n=size(model.node,1);
      m=size(model.member,1);
      r=sum(sum(model.react));
      statdet=false;
      if (m+r==2*n)
        statdet=true;
      end
      if not(statdet)
        warning('The problem is statically indeterminate!')
      end
    end
    
    %%
    function [force,react]=computeMJoints(model)
    % computeMJoints
    % Compute reaction and member forces using method of joints 
      n=size(model.node,1);
      m=size(model.member,1);
      theta=atan2(model.node(model.member(:,2),2)-model.node(model.member(:,1),2),model.node(model.member(:,2),1)-model.node(model.member(:,1),1));
      A=zeros(2*n);
      for i=1:n
        ix=i;
        iy=i+n;
        for j=1:m
          if (model.member(j,1)==i)
            A(ix,j)=cos(theta(j));
            A(iy,j)=sin(theta(j));
          elseif (model.member(j,2)==i)
            A(ix,j)=-cos(theta(j));
            A(iy,j)=-sin(theta(j));
          end
        end
      end
      k=m;
      for i=1:n
        ix=i;
        iy=i+n;
        if model.react(i,1)
          k=k+1;
          A(ix,k)=1;
        end
        if model.react(i,2)
          k=k+1;
          A(iy,k)=1;
        end
      end
      N=-[model.load(:,1); model.load(:,2)];
      Ainv=pinv(A);
      force=Ainv*N;
      [i,j]=find(model.react);
      react=round(sparse(i,j,force(m+1:end),n,2),9);
      force=round(force(1:m),9);
    end
    
    %%
    function delta=computeMVirtualWork(model,force)
    % computeMVirtualWork.m
    % Compute dispacement of nodes using method of virtual work
    % Must precompute forces in members first
      m=size(model.member,1);
      n=size(model.node,1);
      delta=zeros(n,2);
      length=sqrt((model.node(model.member(:,2),2)-model.node(model.member(:,1),2)).^2+(model.node(model.member(:,2),1)-model.node(model.member(:,1),1)).^2);
      for i=1:2*n
        if model.virt(i)
          model.load(:)=0;
          model.load(i)=1;
          virtual=trussJVWFuncs.computeMJoints(model);
          delta(i)=sum(virtual.*force.*length./reshape(model.A.*model.E,m,1));
        end
      end
    end
    
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
      w=ceil(model.A(force>0)./max(model.A)*10);
      for i=1:numel(h)
        set(h(i),'LineWidth',w(i));
      end
      h=plot(reshape(model.node(model.member(force<0,:),1),numel(find(force<0)),2)',reshape(model.node(model.member(force<0,:),2),numel(find(force<0)),2)','b');
      w=ceil(model.A(force<0)./max(model.A)*10);
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
    
  end
end

