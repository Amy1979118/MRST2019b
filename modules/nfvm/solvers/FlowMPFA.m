function [state,matrixInfo] = FlowMPFA(G,T,fluid,W)
%Incompressible single-phase flow using MPFA (MPFA-O, EMPFA, eMPFA)
%  Applicable to general 2D and 3D grids
%  Grid - Grid structure of MRST
%  T - face transmissiblity cell array, one array for each face
%  src - Source structure of MRST
%% Assemble Ax=b
nc=G.cells.num;
nf=G.faces.num;
ncf=100;
b=zeros(nc,1);
[~,rho]=fluid.properties();
[I,J,V]=deal(zeros(ncf*nc,1));k=1;
for i_face=1:nf
    t=T{i_face};
    c1=G.faces.neighbors(i_face,1);
    c2=G.faces.neighbors(i_face,2);
    if(all([c1 c2]~=0))
        ind=t(:,1)~=0;num=sum(ind);
        I(k:k+num-1)=c1;J(k:k+num-1)=t(ind,1);
        V(k:k+num-1)=t(ind,2);k=k+num;
        I(k:k+num-1)=c2;J(k:k+num-1)=t(ind,1);
        V(k:k+num-1)=-t(ind,2);k=k+num;
        if(~all(ind)),b(c1)=b(c1)-t(end,2);b(c2)=b(c2)+t(end,2);end
    else
        c=max(c1,c2);if(c==c2),t(:,2)=-t(:,2);end
        num=size(t,1)-1;
        I(k:k+num-1)=c;J(k:k+num-1)=t(1:end-1,1);
        V(k:k+num-1)=t(1:end-1,2);k=k+num;
         b(c)=b(c)-t(end,2);  % there will always be a constant for boundary faces      
    end
end

%---------------------------------------------------------------------
for i=1:numel(W)
    if(strcmpi(W(i).type,'bhp'))
        pbh=W(i).val;dZ=W(i).dZ;
        for j=1:numel(W(i).cells)
            mycell=W(i).cells(j);
            I(k)=mycell;J(k)=mycell;V(k)=W(i).WI(j);k=k+1;
            b(mycell)=b(mycell)+W(i).WI(j)*(pbh+rho*9.81*dZ(j));
        end 
    else
        myrate=W(i).val.*W(i).WI./sum(W(i).WI);
        mycell=W(i).cells;
        b(mycell)=b(mycell)+myrate;
    end
end
%-------------------------------------------------------------------------
I(k:end)=[];J(k:end)=[];V(k:end)=[];
A=sparse(I,J,V,nc,nc);
if(condest(A)>1e10)
    A(1)=2*A(1);
end
u=A\b;
if(nargout>1)
    matrixInfo.A=A;
    matrixInfo.b=b;
end
%% Compute flux of faces
flux=zeros(nf,1);
for i_face=1:nf
    t=T{i_face};
    if(all(t(:,1)~=0))
        flux(i_face)=u(t(:,1))'*t(:,2);
    else
        flux(i_face)=u(t(1:end-1,1))'*t(1:end-1,2)+t(end,2);
    end
end

state.pressure=u;
state.flux=flux;

%% compute well solutions
wellsol=repmat(struct('pressure',[],'flux',[]),[numel(W) 1]);
for i=1:numel(W)
    if(strcmpi(W(i).type,'bhp'))
        pbh=W(i).val;dZ=W(i).dZ;
        wellsol(i).pressure=pbh+rho*9.81*dZ;
        wellsol(i).flux=W(i).WI.*(wellsol(i).pressure-u(W(i).cells));
    else
        warning('code under development!');
        % write code here babbabababaababababababababababababababababababab
    end
end
state.wellSol=wellsol;
end

