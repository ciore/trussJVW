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
% computeMJoints
% Compute reaction and member forces using method of joints 

function [force,react]=computeMJoints(model)
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