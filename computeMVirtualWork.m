%Compute dispacement of nodes using method of virtual work
%Must precompute forces in members first

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
