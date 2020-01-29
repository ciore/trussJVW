%Compute reaction and member forces using method of joints 

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