%Checks if Nmembers+Nreactions==2*Nnodes

function statdet=checkStaticallyDeterminate(model)
  n=size(model.node,1);
  m=size(model.member,1);
  r=sum(sum(model.react));
  statdet=false;
  if (m+r==2*n)
    statdet=true;
  end
end
