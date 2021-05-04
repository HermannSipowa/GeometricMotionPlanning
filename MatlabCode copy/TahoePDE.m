function [c,f,s] = TahoePDE(x,t,u,DuDx)
% Define PDE; Evaluate right-hand-side of GHF
global Target Chasser Penalty Progress

if t > Progress
    Progress = t
end

n = size(u,1);

%% Defining the Obstacle(s)
lambda  = 1;
plambda = zeros(1,n);


%% Definint the Riemannian Metric
Iinv = inv(diag(Chasser.Moment_Of_Inertia_Calculator()));
F(:,1)    = [u(1:4);zeros(9,1)];
F(:,2:4)  = [zeros(1,3); eye(3); zeros(9,3)];
F(:,5:7)  = [zeros(4,3);Iinv;zeros(6,3)];
F(:,8:13) = [zeros(7,6); eye(6)];
D = diag([Penalty*ones(4,1);ones(3,1);100*Penalty*ones(3,1);10*Penalty*ones(3,1)]);
Finv = inv(F);
H = Finv'*D*Finv; % H = (F')\D/F;
G = lambda*H;


%% Partial derivatives of G
pH_pq0 = -2*Finv'*D*Finv*[[1; zeros(12,1)] zeros(13,12)]*Finv;
pH_pq1 = -2*Finv'*D*Finv*[[0; 1; zeros(11,1)] zeros(13,12)]*Finv;
pH_pq2 = -2*Finv'*D*Finv*[[0; 0; 1; zeros(10,1)] zeros(13,12)]*Finv;
pH_pq3 = -2*Finv'*D*Finv*[[0; 0; 0; 1; zeros(9,1)] zeros(13,12)]*Finv;
for i = 1:n
    pG(:,:,i) = plambda(i)*H;
end
pG(:,:,1) = pG(:,:,1) + lambda*pH_pq0; % pG_pq0
pG(:,:,2) = pG(:,:,2) + lambda*pH_pq1; % pG_pq1
pG(:,:,3) = pG(:,:,3) + lambda*pH_pq2; % pG_pq2
pG(:,:,4) = pG(:,:,4) + lambda*pH_pq3; % pG_pq3


%% Evaluate the christoffels symbols
invG  = inv(G);
Chris = zeros(size(pG));
for i=1:n
    for j=1:n
        for k=1:n
            for l=1:n
                Chris(i,j,k) = Chris(i,j,k) + 1/2*invG(i,l)*(pG(l,j,k)...
                    + pG(l,k,j)-pG(j,k,l));
            end
        end
    end
end


%% fD and its partal derivative
Fd  = RelativeMotionODE(x,u,Target,Chasser);
pFd = Partial_Fd_Quaternion(x,u,Target,Chasser);


%% Computing the coeficients to the PDE
s       = zeros(n,1);
MatProd = zeros(n,1);
for i = 1:n
    MatProd(i) = (DuDx-Fd)'*pG(:,:,i)*Fd;
end
rvec = invG*(pFd'*G*(DuDx-Fd) + 1/2*MatProd);
c = ones(n,1);
f = DuDx-Fd;
for i = 1:n
    s(i) = (DuDx-Fd)'*squeeze(Chris(i,:,:))*(DuDx-Fd)+rvec(i);
end

end
% --------------------------------------------------------------------------