function [p, dpdx, Crtl, sol, t, x, xnew, cost, p_initial] = AGHF(tmax, xpoints, xpoints_new, tpoints, k, X, Xf, T, N, M, AgentNum)
% unicycle_drift_AGHF Solve the AGHF for dynamical unicycle
% code auther: Yinai Fan
%
% system: dot_x = f(x) + F(x)*u
% Lagragian of actuated curve length: (x_dot-f(x))'*G(x_dot-f(x))
%
% The AGHF is equivalent as the Euler Lagrange equation:
% dx/ds = AGHF  <===> dx/ds = - (dL/dx - d(dL/d(dot_x))/dt)
% In this implementation, the EL equation is used for solving PDE
%
% Inputs:
% tmax -- smax
% xpoints -- number of grids in time t
% intgrids -- number of grids for integration of solution
% tpoints -- number of grids in s
% k -- penalty (lambda is used in paper)
% X0, Xf -- initial and final value of states
% T -- motion duration
%
% Outputs:
% p, dp -- steady state PDE solution for x, dot_x of the AGHF
% xfull -- state trajectory integrated by using full control (virtual control included)
% xdrift -- solution of the motion, state trajectory integrated by using actual control
% bufull -- full control (virtual control included)
% budrift -- actual control
% sol -- raw solution of pdepe
% t -- s grids vector for solving AGHF,  x -- time grids vector for solving AGHF
% xnew -- time vector for integration of solutions
% cost -- actuated curve length values along s grids

%----------------------- pdepe setup ------------------------------
m = 0;
x = linspace(0,T,xpoints);                  % discretization of the curve in time
t = [0 logspace(-4,log10(tmax),tpoints-1)]; % discretization of s interval in log scale

% global variable for the Euler Lagrange function and Metric G
global get_EL get_G Period

% generate functions for EL, metric G, drift term f, control matrix F
% please modify this function if system is changed
[get_EL, get_G, get_f, get_F] = generate_fh_of_model(N,M,AgentNum); 

% opts = odeset('RelTol',2.22045e-9,'AbsTol',2.22045e-20);
opts = odeset('AbsTol',1e-14);
tic;
disp('solving PDE...');
% Solve AGHF
sol = pdepe(m,@(x,t,u,DuDx) mypdexpde(x,t,u,DuDx,k,N),...
              @(x) mypdexic(x, X, Xf, T),...
              @(xl,ul,xr,ur,t) mypdexbc(xl,ul,xr,ur,t, X, Xf, N),...
              x,t,opts);
% The solution "sol" is of form sol(t,x,i)
toc;

%------------------------ Control Extraction ------------------------------
tic;
disp('Extracting Trajectory...');

% initialize controls and states
Crtl = struct('u', cell(1, AgentNum)); % actual control
xnew = linspace(0,T,xpoints_new); % time grids for integration

% use spline to interpolate state from steady state pde solution, p
p         = zeros(N,xpoints_new);
p_initial = zeros(N,xpoints_new);
dpdx      = zeros(N,xpoints_new);
dx_new    = diff([xnew(3),xnew,xnew(xpoints_new-2)]);
for i = 1:N
    p(i,:)         = spline(x,sol(end,:,i),xnew);
    p_initial(i,:) = spline(x,sol(1,:,i),xnew);
    
    % Numerically computing the time derivative of the state
    dp_new    = diff([p(i,3),p(i,:),p(i,xpoints_new-2)]);
    dpdx(i,:) = (dp_new(1:end-1)./dx_new(1:end-1).*dx_new(2:end) ...
        + dp_new(2:end)./dx_new(2:end).*dx_new(1:end-1)) ...
        ./ (dx_new(2:end)+dx_new(1:end-1));
end

% control extraction
mm = 6;
for i = 1 : xpoints_new
    % f(x(t)), note: the F_d(x) in paper
    drift = get_f(p(:,i));
    % get [Fc F], note: the F_bar matrix in paper
    Ffull = get_F(p(:,i));
    B = Ffull(:,M+1:end);
    
    % Extracting/approximate the required control
    for jj = 1:AgentNum
        idx = 1+mm*(jj-1):mm*jj;
        Crtl(jj).u(:,i) = (B.'*B)\(B.')* ( dpdx(idx,i) - drift(idx) );
    end
end

% Redimensionalizing the problem
p_initial = Period.*p_initial;
p         = Period.*p;


%------------------ Calculating the action integral -----------------------
% initialize cost -- actuated curve length
cost = zeros(1,tpoints);
% calculate the atuated curve length for each grid of s
X_temp = zeros(N,xpoints);
dX_temp = zeros(N,xpoints);
for j = 1:tpoints
    for kk = 1:N
    [X_temp(kk,:), dX_temp(kk,:)]=pdeval(m,x,sol(j,:,kk),x);
    end
    %for each s grid, integrate  (x_dot-f(x))'*G(x_dot-f(x)) to get curve
    %length
    for i = 1:xpoints
        X = X_temp(:,i);
        dX = dX_temp(:,i);
        f = get_f(X);
        G = get_G(X,k);
        cost(j) = cost(j) + (dX-f)'*G*(dX-f)*(T/(xpoints-1));
    end
end
toc;
disp('Done!!!!!');

end
% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------

% The followings are the PDE, initial condition and boundary condition for pdepe:

function [c,f,s] = mypdexpde(x,t,u,DuDx,k,N)    % Define PDE; right-hand-side of AGHF

global get_EL

% evaluate EL
EL = get_EL(u,DuDx,k);
% dL/dx
pLx = EL(:,1);
% dL/d(dot_x)
pLxd = EL(:,2);

f = pLxd;
s = -pLx;

c = ones(N,1);

% for details of f, s, c, please see the documentation for pdepe

end
% -------------------------------------------------------------------------
 
function u0 = mypdexic(x, X0, Xf, T)    % Initial condition for AGHF

% % straight line connecting X0 and Xf
% u0=X0+(Xf-X0)*(x/T);
% % add some sinusoidal deviation to x1 initial guess (can be discarded)
% u0(1)=X0(1)+(Xf(1)-X0(1))*(x/T) + 0.01*sin(2*pi*x/T);

%% Sinusoidal initial condition
freq1 = pi/2 + 2*pi;
freq2 = pi + 2*pi;

u0 = X0*cos(freq1*x/T) ...
    + Xf*sin(freq1*x/T) ...
    - (X0+Xf)*sin(freq2*x/T);

end
% -------------------------------------------------------------------------

function [pl,ql,pr,qr] = mypdexbc(xl,ul,xr,ur,t, X0, Xf, N)    % Boundary condition for AGHF

pl = ul-X0;
ql = zeros(N,1);
pr = ur-Xf;
qr = zeros(N,1);

end
% -------------------------------------------------------------------------

% generate functions for calculating EL equation:
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!! Need to modify this function if system is changed !!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Inputs:
% N,M -- # of states and inputs
% Output:
% get_EL -- function handle for calculationg EL equation, it returns dL/dx
% and % dL/d(dot_x), it takes x and x_dot as input
% get_G -- function handle for calculating the metric G(x)
% get_f -- function handle for calculating the drift term f(x)
% get_F -- function handle for calculating the full control matrix [Fc(x) F(x)]

function [get_EL, get_G, get_f, get_F] = generate_fh_of_model(N,M,AgentNum) 
tic;
disp('generating functions for pdepe...');

% ------------------------- system setup ----------------------------------
syms Penalty
global mu_dot Period
A  = [zeros(3) eye(3);
      3*mu_dot^2 0 0 0 2*mu_dot 0;
      0 0 0 -2*mu_dot 0 0;
      0 0 -mu_dot^2 0 0 0]*Period; % Partial derivative pFd the drift vector field
N = 6;
  % make the state and state deirvative simbolic
  for i = 1:AgentNum
      X_Variables{i} = [{strcat('x',num2str(i))},{strcat('y',num2str(i))},...
          {strcat('z',num2str(i))},{strcat('xdot',num2str(i))},...
          {strcat('ydot',num2str(i))},{strcat('zdot',num2str(i))}];
      dX_Variables{i} = [{strcat('dx',num2str(i))},{strcat('dy',num2str(i))},...
          {strcat('dz',num2str(i))},{strcat('dxdot',num2str(i))},...
          {strcat('dydot',num2str(i))},{strcat('dzdot',num2str(i))}];
      if i == 1
          AMat_Aug = A;
          % [Fc F], note: the F_bar matrix in paper
          F = sym(eye(N)); F_Aug = F;
          % penalty matrix (D matrix in the paper)
          D = diag([Penalty*ones(1,N-M) ones(1,M)]); D_Aug = D;
      else
          AMat_Aug = blkdiag(AMat_Aug,A);
          F_Aug = blkdiag(F_Aug,F);
          D_Aug = blkdiag(D_Aug,D);
      end
  end
% syms x y z xdot ydot zdot  dx dy dz dxdot dydot dzdot Penalty real
X1  = cell2sym( X_Variables{1} );
X2  = cell2sym( X_Variables{2} );
X3  = cell2sym( X_Variables{3} );
X4  = cell2sym( X_Variables{4} );
% X5  = cell2sym( X_Variables{5} );
% X6  = cell2sym( X_Variables{6} );

dX1 = cell2sym( dX_Variables{1} );
dX2 = cell2sym( dX_Variables{1} );
dX3 = cell2sym( dX_Variables{1} );
dX4 = cell2sym( dX_Variables{1} );
% dX5 = cell2sym( dX_Variables{1} );
% dX6 = cell2sym( dX_Variables{1} );

Xaug  = [X1(:);X2(:);X3(:);X4(:)]; % ;X5(:);X6(:)];
dXaug = [dX1(:);dX2(:);dX3(:);dX4(:)]; % ;dX5(:);dX6(:)];


% drift term f(x) note: Fd(x) is used in the paper
f = AMat_Aug*Xaug;
% the metric, before adding the state constraint barrier functions
H = (F_Aug.')^(-1)*D_Aug*F_Aug^(-1);


% -------------------- state constraint barrier function ------------------
% B is a vector of barrier function for 1 state constraint, each scaler state
% constraint need 1 barrier function
B = [];
% barrier function parameters
% kb = .01; % barrier function gain
% pb = 1; % order of barrier function

% b is one penalty term for state constraint
b = 0; % no state contraints
% b = kb/((2*pi/4)^2 - z2^2)^pb; %contraint on angular velocity, |z2|<pi/2
% b = kb/(2^2 - z1^2)^pb; %contraint on translational velocity, |z1|<2
% b = kb/(z1+2)^pb + kb/(2-z1)^pb; %contraint on translational velocity |z1|<2

B = [B b]; % attach the barrier function to B

% -------------------- Metric and curve Length ----------------------------

% the metric with state constraint barier functions
G = (sum(B)+1)*H;
% actuated curve length
L = simplify( (dXaug - f).' * G * (dXaug - f) );

% -------------------- Function Generations -------------------------------

% taking derivatives symbolically to get the EL euqation terms
% pLx -- dL/dx
% pLxd -- dL/d(x_dot)
pLx  = sym('pLx', [N*AgentNum 1],'real');
pLxd = sym('pLxd', [N*AgentNum 1],'real');
for i=1:N*AgentNum
    pLx(i)  = diff(L,Xaug(i));
    pLxd(i) = diff(L,dXaug(i));
end

% generate functions
get_EL = matlabFunction([pLx, pLxd],'vars', {Xaug,dXaug,Penalty});
get_G  = matlabFunction(G,'vars', {Xaug,Penalty});
get_f  = matlabFunction(f,'vars',{Xaug});
get_F  = matlabFunction(F,'vars',{Xaug});

toc;
end


function [xdot] = MyOde(t,x,u)
format long 
global mu_dot Period
% Reshaping the input into an 6 by n matrix
[l,~] = size(x);
x = reshape(x,6,l/6);

% Computing the CW control and dynamic matrices
B = [zeros(3); eye(3)];
A = [zeros(3) eye(3);
    3*mu_dot^2 0 0 0 2*mu_dot 0;
    0 0 0 -2*mu_dot 0 0;
    0 0 -mu_dot^2 0 0 0]*Period;

% Computing the derivatives for the ode
xdot = A*x + B*u;
xdot = xdot(:);

end
% -------------------------------------------------------------------------






