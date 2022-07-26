clear all
close all
rng(2141444)


%% *************************** Dynamics ***********************************

f_u =  @ChemostatMonodKinetcs;
n = 2;
m = 1; % number of control inputs


%% ************************** Discretization ******************************

deltaT = 0.1;
%Runge-Kutta 4
k1 = @(t,x,u) (  f_u(t,x,u) );
k2 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT/2,u) );
k3 = @(t,x,u) ( f_u(t,x + k2(t,x,u)*deltaT/2,u) );
k4 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT,u) );
f_ud = @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   );

%% ************************** Basis functions *****************************

basisFunction = 'rbf';
% RBF centers
Nrbf = 20;
cent = rand(n,Nrbf)*20;
rbf_type = 'thinplate' ; %'polyharmonic' 'invmultquad' 'invquad' 'gauss' 'thinplate'
% Lifting mapping - RBFs + the state itself
liftFun = @(xx)( [xx;rbf(xx,cent,rbf_type)] );
Nlift = Nrbf + n;


%% ************************** Collect data ********************************
tic
disp('Starting data collection')
Nsim = 1000;
Ntraj = 200;
Ubig = rand([Nsim Ntraj])/2;
 
% Random initial conditions
X1current = 5*((rand(1,Ntraj)*2)+2);
X2current = 2*((rand(1,Ntraj)*4)+1);
Xcurrent = [X1current;X2current];

X = []; Y = []; U = [];
for i = 1:Nsim
    Xnext = f_ud(0,Xcurrent,Ubig(i,:));
    X = [X Xcurrent];
    Y = [Y Xnext];
    U = [U Ubig(i,:)];
    Xcurrent = Xnext;
end
fprintf('Data collection DONE, time = %1.2f s \n', toc);


%% ******************************* Lift ***********************************

disp('Starting LIFTING')
tic
Xlift = liftFun(X);
Ylift = liftFun(Y);
fprintf('Lifting DONE, time = %1.2f s \n', toc);

%% ********************** Build predictor *********************************

disp('Starting REGRESSION')
tic
W = [Ylift ; X];
V = [Xlift; U];
VVt = V*V';
WVt = W*V';
M = WVt * pinv(VVt); % Matrix [A B; C 0]
Alift = M(1:Nlift,1:Nlift);
Blift = M(1:Nlift,Nlift+1:end);
Clift = M(Nlift+1:end,1:Nlift);

fprintf('Regression done, time = %1.2f s \n', toc);

%% *********************** Predictor comparison ***************************

Tmax = 30;
Nsim = Tmax/deltaT;
uprbs = myprbs(Nsim,.5)/4+.1;
u_dt = @(i)(  uprbs(i+1) );
f_cont_d = @(t,xx)( f_ud(t,xx,u_dt(t)) );

%Initial condition
%X1 =5*((rand(1)*2)+2);
%X2 =2*((rand(1)*4)+1);
xe = [4.1459;15.8541]; % This is the equlibrium of the chemostat model with the used parameters
X1 = 15;
X2 = 7;
x0 = [X1;X2];
x_true = x0;
% Lifted initial condition
xlift = liftFun(x0);

% Local linearization predictor at x0
x = sym('x',[2;1]); u = sym('u',[1;1]);
Ac_x0 = double(subs(jacobian(f_u(0,x,u),x),[x;u],[x0;0]));
Bc_x0 = double(subs(jacobian(f_u(0,x,u),u),[x;u],[x0;0]));
c_x0 = double(subs(f_u(0,x,u),[x;u],[x0;0])) - Ac_x0*x0 - Bc_x0*0;
ABc = expm([Ac_x0 [Bc_x0 c_x0] ; zeros(2,4)]*deltaT); % discretize
Ad_x0 = ABc(1:2,1:2); Bd_x0 = ABc(1:2,3); cd_x0 = ABc(1:2,4);
X_loc_x0 = x0;

% Local linearization predictor at xe
x = sym('x',[2;1]); u = sym('u',[1;1]);
Ac_xe = double(subs(jacobian(f_u(0,x,u),x),[x;u],[xe;0]));
Bc_xe = double(subs(jacobian(f_u(0,x,u),u),[x;u],[xe;0]));
c_xe = double(subs(f_u(0,x,u),[x;u],[xe;0])) - Ac_xe*xe - Bc_xe*0;
ABc = expm([Ac_xe [Bc_xe c_xe] ; zeros(2,4)]*deltaT); % discretize
Ad_xe = ABc(1:2,1:2); Bd_xe = ABc(1:2,3); cd_xe = ABc(1:2,4); 
X_loc_xe = xe;


% Simulate
for i = 0:Nsim-1
    % Koopman predictor
    xlift = [xlift, Alift*xlift(:,end) + Blift*u_dt(i)]; % Lifted dynamics
    
    % True dynamics
    x_true = [x_true, f_ud(0,x_true(:,end),u_dt(i)) ];
    % Local linearization predictor at x0
    X_loc_x0 = [X_loc_x0, Ad_x0*X_loc_x0(:,end) + Bd_x0*u_dt(i) + cd_x0];
    
    % Local linearization predictor at xe
    X_loc_xe = [X_loc_xe, Ad_xe*X_loc_xe(:,end) + Bd_xe*u_dt(i) + cd_xe];
    
end
x_koop = Clift * xlift; % Koopman predictions

fprintf('Local at x0 Err %f \n',  rmse(x_true,X_loc_x0));
fprintf('Local at xe Err %f \n',  rmse(x_true,X_loc_xe));
fprintf('KoopmanErr %f \n',  rmse(x_true,x_koop)) ;

%% ****************************  Plots  ***********************************

lw = 4;

figure
stairs((0:Nsim-1)*deltaT,u_dt(0:Nsim-1),'linewidth',lw); hold on
title('Control input $u = D$', 'interpreter','latex'); xlabel('Time [d]','interpreter','latex');
set(gca,'fontsize',20)

figure
plot((0:Nsim)*deltaT,x_true(2,:),'linewidth',lw); hold on
plot((0:Nsim)*deltaT,x_koop(2,:), '--r','linewidth',lw)
plot((0:Nsim)*deltaT,X_loc_x0(2,:), '--g','linewidth',lw-1)
plot((0:Nsim)*deltaT,X_loc_xe(2,:), '--k','linewidth',lw-1)
axis([0 Tmax min(x_koop(2,:))-0.15 max(x_koop(2,:))+0.15])
title('Model comparison for biomass concentration $X$','interpreter','latex'); xlabel('Time [d]','interpreter','latex');
set(gca,'fontsize',20)
LEG = legend('Real system','Koopman model','Local at $x_0$','Local at $xe$','location','southwest');
set(LEG,'interpreter','latex')

figure
plot((0:Nsim)*deltaT,x_true(1,:),'linewidth',lw); hold on
plot((0:Nsim)*deltaT,x_koop(1,:), '--r','linewidth',lw)
plot((0:Nsim)*deltaT,X_loc_x0(1,:), '--g','linewidth',lw-1)
plot((0:Nsim)*deltaT,X_loc_xe(1,:), '--k','linewidth',lw-1)
axis([0 Tmax min(x_koop(1,:))-0.1 max(x_koop(1,:))+0.1])
title('Model comparison for   substrate concentration $S$','interpreter','latex'); xlabel('Time [d]','interpreter','latex');
set(gca,'fontsize',20)
LEG = legend('Real system','Koopman model','Local at $x_0$','Local at $xe$','location','southwest');
set(LEG,'interpreter','latex')

%%
% RMSE comparison 
disp('Model Compartision with RMSE')
fprintf('x10 = S0 %f \n',  X1);
fprintf('x20 = X0 %f \n',  X2);
fprintf('Local at x0 Err %f \n',  rmse(x_true,X_loc_x0));
fprintf('Local at xe Err %f \n',  rmse(x_true,X_loc_xe));
fprintf('KoopmanErr %f \n',  rmse(x_true,x_koop)) ;
