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
% Checks if Nmembers+Nreactions==2*Nnodes

function statdet=checkStaticallyDeterminate(model)
  n=size(model.node,1);
  m=size(model.member,1);
  r=sum(sum(model.react));
  statdet=false;
  if (m+r==2*n)
    statdet=true;
  end
end
