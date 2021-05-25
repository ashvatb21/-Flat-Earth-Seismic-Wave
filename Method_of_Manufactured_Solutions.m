%% Spatial Convergence MMS Trapezoidal

clear all, close all, clc

%params for problem
a = 0; % Lower Bound x (km)
b = 1; % Upper Bound x (km)
w = 3.25; % Frequency
k = 7.5; % Wave Number
c = w/k; % Wave Speed (km/s)
Amp = 15; % Amplitude (km)
ln = b-a; % The domian of x (km)
T = 1; % Final time to run to (seconds)
dt = 1e-3; % Make dt small in spatial convergence test so that time error doesn't pollute convergence

%Initial Conditions, Gaussian Function
uex = @(x,t) Amp*sin(w*t - k*x); % Exact Solution
u_init = @(x) uex(x,0); % Initial Displacement
v_init = @(x) w*Amp*cos(-k*x); % Initial Velocity

%g(x,t) and Boundary Conditions
bc_a = @(t) Amp*sin(w*t - k*a); % Dirichlet Boundary Condition at 'a'
bc_b = @(t) Amp*sin(w*t - k*b); % Dirichlet Boundary Condition at 'b'
fcn = @(x,t) c^2.*Amp.*k^2.*sin(w*t-k*x) - Amp.*w^2.*sin(w*t-k*x); % Forcing Function

%# of n points to use
nvect = [40; 80; 120; 240];

%initialize error
err = zeros( size( nvect ) );

tic
for j = 1 : length( nvect )
    
    n = nvect(j);
    
    %---Build n, xj points, A matrix and g vector
    
    %Build interp points
    xj = ( a + (b-a)*(0:n)./n )';
    
    %grid spacing (uniform)
    dx = (b-a)/n;
    
    %Build A matrix
    %Use truncated version from lecture notes
    A_temp = (c^2/(dx^2))*(-2*diag( ones(n-1,1),0) + 1*diag( ones(n-2,1),-1) + ...
        1*diag( ones(n-2,1),1) );
    
    % Matrix of Zeros of size A_temp
    Zeros = zeros(size(A_temp));
    
    % The expanded A Matrix for 2nd degree differential equation
    A = [Zeros eye(size(A_temp)); A_temp Zeros];
    
    % Identity matrix of size A
    I = eye(size(A));
    
    %Build g for this set of xj
    g1 = @(t) [zeros(size( xj(2:end-1) ));fcn(xj(2:end-1),t)];
    g2 = @(t) [zeros(size(xj(2:end-1))); c^2/dx^2*bc_a(t); zeros(size(xj(3:end-2))); c^2/dx^2*bc_b(t)];
    g = @(t) g1(t) + g2(t);
    
    %Build RHS for IVP, f(u,t)
    f = @(u,t) A*u + g(t);
    %---
    
    %---Initialize for time stepping
    uk = [u_init(xj(2:end-1));v_init(xj(2:end-1))];
    tk = 0;
    tvect = dt : dt : T;
    
    %# snapshots to save
    nsnps = 250;
    ind = max( 1, round(length(tvect)/nsnps) );
    tsv = tvect( 1 : ind : end );
    
    u = zeros( 2*n-2, length(tsv));
    
    cnt = 1;
    %---
    
    %---Do time stepping
    for jj = 1 : length( tvect )
        
        tkp1 = tk + dt;
        
        %Update solution at next time using trap method
        ukp1 = (((I-dt/2*A)) \ (uk + dt/2*( f(uk,tk)  + g(tkp1) )));
        
        uk = ukp1;
        tk = tkp1;
        
        if min(abs( tkp1-tsv ) ) < 1e-8
            
            u(:,cnt) = uk;
            cnt = cnt + 1;
            
        end
        
    end
    err(j) = norm( ukp1(1:n-1) - uex(xj(2:end-1),T) ) / norm( uex(xj(2:end-1),T) );
end
toc

[X,T] = meshgrid( xj(2:end-1), tsv );
U = u(1:n-1,:).';

% -Waterfall plot of solns

figure(1)
subplot(1,2,1)
waterfall( X,T,U )
set( gca, 'fontsize', 15, 'ticklabelinterpreter', 'latex' )
title('$u_{FD}$', 'fontsize', 20, 'interpreter', 'latex')
xlabel('$x$', 'fontsize', 15, 'interpreter', 'latex')
ylabel('$t$', 'fontsize', 15, 'interpreter', 'latex')

subplot(1,2,2)
waterfall( X,T, uex(X,T) )
set( gca, 'fontsize', 15, 'ticklabelinterpreter', 'latex' )
title('$u_{exact}$', 'fontsize', 20, 'interpreter', 'latex')
xlabel('$x$', 'fontsize', 15, 'interpreter', 'latex')
ylabel('$t$', 'fontsize', 15, 'interpreter', 'latex')
%

%--Error plots
figure(2)
c1 = err(end)/(dx^2);
loglog( ln./nvect, c1*(ln./nvect).^2, 'k--', 'linewidth', 2 ), hold on
title('Spatial Convergence Trapezoidal')

%plot err
loglog( ln./nvect, err , 'b.', 'markersize', 26 )
h = legend('$O(\Delta x^2)$', '$Error$');
set(h, 'Interpreter','latex', 'fontsize', 16, 'Location', 'NorthWest' )
%make pretty
xlabel( '$\Delta x$', 'interpreter', 'latex', 'fontsize', 16)
ylabel( '$||\textbf{e}||/||\textbf{u}_e||$ ', 'interpreter', 'latex', 'fontsize', 16)
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16 )


%% Temporial Convergence MMS Trapezoidal

clear all, close all, clc

%params for problem 
a = 0; % Lower Bound x (km)
b = 1; % Upper Bound x (km)
w = 3.25; % Frequency
k = 7.5; % Wave Number
c = w/k; % Wave Speed (km/s)
Amp = 15; % Amplitude (km)
ln = b-a; % The domian of x (km)
T = 1; % Final time to run to (seconds)
dt = 1e-3; % Make dt small in spatial convergence test so that time error doesn't pollute convergence
n = 1000;

%Initial Conditions, Gaussian Function
uex = @(x,t) Amp*sin(w*t - k*x); % Exact Solution
u_init = @(x) uex(x,0); % Initial Displacement
v_init = @(x) w*Amp*cos(-k*x); % Initial Velocity

%g(x,t) and Boundary Conditions
bc_a = @(t) Amp*sin(w*t - k*a); % Dirichlet Boundary Condition at 'a'
bc_b = @(t) Amp*sin(w*t - k*b); % Dirichlet Boundary Condition at 'b'
fcn = @(x,t) c^2.*Amp.*k^2.*sin(w*t-k*x) - Amp.*w^2.*sin(w*t-k*x); % Forcing Function

%# of n points to use
dtvect = [2.5; 1; 5e-1; 2.5e-1; 1e-1]*0.25; % From HW8

%initialize error
err = zeros(size(dtvect));

tic
for j = 1:length(dtvect)
    
    dt = dtvect(j);
    
    %Build interp points
    xj = ( a + (b-a)*(0:n)./n )';
    
    %grid spacing (uniform)
    dx = (b-a)/n;
    
    %Build A matrix
    %Use truncated version from lecture notes
    A_temp = (c^2/(dx^2))*(-2*diag( ones(n-1,1),0) + 1*diag( ones(n-2,1),-1) + ...
        1*diag( ones(n-2,1),1) );
    
    % Matrix of Zeros of size A_temp
    Zeros = zeros(size(A_temp));
    
    % The expanded A Matrix for 2nd degree differential equation
    A = [Zeros eye(size(A_temp)); A_temp Zeros];
    
    % Identity matrix of size A
    I = eye(size(A));
    
    %Build g for this set of xj
    g1 = @(t) [zeros(size( xj(2:end-1) ));fcn(xj(2:end-1),t)];
    g2 = @(t) [zeros(size(xj(2:end-1))); c^2/dx^2*bc_a(t); zeros(size(xj(3:end-2))); c^2/dx^2*bc_b(t)];
    g = @(t) g1(t) + g2(t);
    
    %Build RHS for IVP, f(u,t)
    f = @(u,t) A*u + g(t);
    %---
    
    %---Initialize for time stepping
    uk = [u_init(xj(2:end-1));v_init(xj(2:end-1))];
    tk = 0;
    tvect = dt : dt : T;
    
    %# snapshots to save
    nsnps = 250;
    ind = max( 1, round(length(tvect)/nsnps) );
    tsv = tvect( 1 : ind : end );
    
    u = zeros( 2*n-2, length(tsv));
    
    cnt = 1;
    %---
    
    %---Do time stepping
    for jj = 1 : length( tvect )
        
        tkp1 = tk + dt;
        
        %Update solution at next time using trap method
        ukp1 = (((I-dt/2*A)) \ (uk + dt/2*( f(uk,tk)  + g(tkp1) )));
        
        uk = ukp1;
        tk = tkp1;
        
        if min(abs( tkp1-tsv ) ) < 1e-8
            
            u(:,cnt) = uk;
            cnt = cnt + 1;
            
        end
        
    end
    err(j) = norm( ukp1(1:n-1) - uex(xj(2:end-1),T) ) / norm( uex(xj(2:end-1),T) );
end
toc

%Error plots
figure(3)
c = err(end)*(1./dt^2);
loglog( dtvect, c*(dtvect).^2, 'k--', 'linewidth', 2 ), hold on

%plot err
loglog(dtvect, err , 'b.', 'markersize', 26 )
title('Temporal Convergence Trapezoidal')

%make pretty
h = legend('$O(\Delta t^2)$', '$Error$');
set(h, 'Interpreter','latex', 'fontsize', 16, 'Location', 'NorthWest' )
xlabel( '$\Delta t$', 'interpreter', 'latex', 'fontsize', 16)
ylabel( '$||\textbf{e}||/||\textbf{u}_e||$ ', 'interpreter', 'latex', 'fontsize', 16)
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16 )

%% Spatial Convergence MMS RK4

clear all, close all, clc

%params for problem
a = 0; % Lower Bound x (km)
b = 1; % Upper Bound x (km)
w = 3.25; % Frequency
k = 7.5; % Wave Number
c = w/k; % Wave Speed (km/s)
Amp = 15; % Amplitude (km)
ln = b-a; % The domian of x (km)
T = 1; % Final time to run to (seconds)
dt = 1e-5; % Make dt small in spatial convergence test so that time error doesn't pollute convergence

%Initial Conditions, Gaussian Function
uex = @(x,t) Amp*sin(w*t - k*x); % Exact Solution
u_init = @(x) uex(x,0); % Initial Displacement
v_init = @(x) w*Amp*cos(-k*x); % Initial Velocity

%g(x,t) and Boundary Conditions
bc_a = @(t) Amp*sin(w*t - k*a); % Dirichlet Boundary Condition at 'a'
bc_b = @(t) Amp*sin(w*t - k*b); % Dirichlet Boundary Condition at 'b'
fcn = @(x,t) c^2.*Amp.*k^2.*sin(w*t-k*x) - Amp.*w^2.*sin(w*t-k*x); % Forcing Function

%# of n points to use
nvect = [40; 80; 120; 240];

%initialize error
err = zeros( size( nvect ) );

tic
for j = 1 : length( nvect )
    
    n = nvect(j);
    
    %---Build n, xj points, A matrix and g vector
    
    %Build interp points
    xj = ( a + (b-a)*(0:n)./n )';
    
    %grid spacing (uniform)
    dx = (b-a)/n;
    
    %Build A matrix
    %Use truncated version from lecture notes
    A_temp = (c^2/(dx^2))*(-2*diag( ones(n-1,1),0) + 1*diag( ones(n-2,1),-1) + ...
        1*diag( ones(n-2,1),1) );
    
    % Matrix of Zeros of size A_temp
    Zeros = zeros(size(A_temp));
    
    % The expanded A Matrix for 2nd degree differential equation
    A = [Zeros eye(size(A_temp)); A_temp Zeros];
    
    % Identity matrix of size A
    I = eye(size(A));
    
    %Build g for this set of xj
    g1 = @(t) [zeros(size( xj(2:end-1) ));fcn(xj(2:end-1),t)];
    g2 = @(t) [zeros(size(xj(2:end-1))); c^2/dx^2*bc_a(t); zeros(size(xj(3:end-2))); c^2/dx^2*bc_b(t)];
    g = @(t) g1(t) + g2(t);
    
    %Build RHS for IVP, f(u,t)
    f = @(u,t) A*u + g(t);
    %---
    
    %---Initialize for time stepping
    uk = [u_init(xj(2:end-1));v_init(xj(2:end-1))];
    tk = 0;
    tvect = dt : dt : T;
    
    %# snapshots to save
    nsnps = 250;
    ind = max( 1, round(length(tvect)/nsnps) );
    tsv = tvect( 1 : ind : end );
    
    u = zeros( n-1, length(tsv));
    
    cnt = 1;
    %---
    
    %---Do time stepping
    for jj = 1:length(tvect)
        
        tkp1 = tk + dt;
        
        y1 = f(uk, tk);
        y2 = f(uk + 0.5*dt*y1, tk);
        y3 = f(uk + 0.5*dt*y2, tk);
        y4 = f(uk + dt*y3, tk);
        
        ukp1 = uk + (1/6)*dt*(y1 + 2*y2 + 2*y3 + y4);
        
        uk = ukp1;
        tk = tkp1;
        
        if min(abs( tkp1-tsv ) ) < 1e-8
            
            u(:,cnt) = uk(1:n-1);
            cnt = cnt + 1;
        end
        
    end
    %---
    err(j) = norm( ukp1(1:n-1) - uex(xj(2:end-1),T) ) / norm( uex(xj(2:end-1),T) );
end
toc

[X,T] = meshgrid( xj(2:end-1), tsv );

% -Waterfall plot of solns

figure(1)
subplot(1,2,1)
waterfall( X,T, u.' )
set( gca, 'fontsize', 15, 'ticklabelinterpreter', 'latex' )
title('$u_{FD}$', 'fontsize', 20, 'interpreter', 'latex')
xlabel('$x$', 'fontsize', 15, 'interpreter', 'latex')
ylabel('$t$', 'fontsize', 15, 'interpreter', 'latex')

subplot(1,2,2)
waterfall( X,T, uex(X,T) )
set( gca, 'fontsize', 15, 'ticklabelinterpreter', 'latex' )
title('$u_{exact}$', 'fontsize', 20, 'interpreter', 'latex')
xlabel('$x$', 'fontsize', 15, 'interpreter', 'latex')
ylabel('$t$', 'fontsize', 15, 'interpreter', 'latex')
%

%--Error plots
figure(2)
c1 = err(end)/(dx^2);
loglog( ln./nvect, c1*(ln./nvect).^2, 'k--', 'linewidth', 2 ), hold on
title('Spatial Convergence RK4')

%plot err
loglog( ln./nvect, err , 'b.', 'markersize', 26 )
h = legend('$O(\Delta x^2)$', '$Error$');
set(h, 'Interpreter','latex', 'fontsize', 16, 'Location', 'NorthWest' )
%make pretty
xlabel( '$\Delta x$', 'interpreter', 'latex', 'fontsize', 16)
ylabel( '$||\textbf{e}||/||\textbf{u}_e||$ ', 'interpreter', 'latex', 'fontsize', 16)
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16 )

%% Temporial Convergence MMS RK4

clear all, close all, clc

%params for problem
a = 0; % Lower Bound x (km)
b = 1; % Upper Bound x (km)
w = 3.25; % Frequency
k = 7.5; % Wave Number
c = w/k; % Wave Speed (km/s)
Amp = 15; % Amplitude (km)
ln = b-a; % The domian of x (km)
T = 1; % Final time to run to (seconds)
dt = 1e-3; % Make dt small in spatial convergence test so that time error doesn't pollute convergence
n = 5000;

%Initial Conditions, Gaussian Function
uex = @(x,t) Amp*sin(w*t - k*x); % Exact Solution
u_init = @(x) uex(x,0); % Initial Displacement
v_init = @(x) w*Amp*cos(-k*x); % Initial Velocity

%g(x,t) and Boundary Conditions
bc_a = @(t) Amp*sin(w*t - k*a); % Dirichlet Boundary Condition at 'a'
bc_b = @(t) Amp*sin(w*t - k*b); % Dirichlet Boundary Condition at 'b'
fcn = @(x,t) c^2.*Amp.*k^2.*sin(w*t-k*x) - Amp.*w^2.*sin(w*t-k*x); % Forcing Function

%# of n points to use
dtvect = [2.5; 1; 5e-1; 2.5e-1; 1e-1]*0.25; % From HW8

%initialize error
err = zeros(size(dtvect));

tic
for j = 1:length(dtvect)
    
    dt = dtvect(j);
    
    %Build interp points
    xj = ( a + (b-a)*(0:n)./n )';
    
    %grid spacing (uniform)
    dx = (b-a)/n;
    
    %Build A matrix
    %Use truncated version from lecture notes
    A_temp = (c^2/(dx^2))*(-2*diag( ones(n-1,1),0) + 1*diag( ones(n-2,1),-1) + ...
        1*diag( ones(n-2,1),1) );
    
    % Matrix of Zeros of size A_temp
    Zeros = zeros(size(A_temp));
    
    % The expanded A Matrix for 2nd degree differential equation
    A = [Zeros eye(size(A_temp)); A_temp Zeros];
    
    % Identity matrix of size A
    I = eye(size(A));
    
    %Build g for this set of xj
    g1 = @(t) [zeros(size( xj(2:end-1) ));fcn(xj(2:end-1),t)];
    g2 = @(t) [zeros(size(xj(2:end-1))); c^2/dx^2*bc_a(t); zeros(size(xj(3:end-2))); c^2/dx^2*bc_b(t)];
    g = @(t) g1(t) + g2(t);
    
    %Build RHS for IVP, f(u,t)
    f = @(u,t) A*u + g(t);
    %---
    
    %---Initialize for time stepping
    uk = [u_init(xj(2:end-1));v_init(xj(2:end-1))];
    tk = 0;
    tvect = dt : dt : T;
    
    %# snapshots to save
    nsnps = 250;
    ind = max( 1, round(length(tvect)/nsnps) );
    tsv = tvect( 1 : ind : end );
    
    u = zeros( n-1, length(tsv));
    
    cnt = 1;
    %---
    
    %---Do time stepping
    for jj = 1:length(tvect)
        
        tkp1 = tk + dt;
        
        y1 = f(uk, tk);
        y2 = f(uk + 0.5*dt*y1, tk);
        y3 = f(uk + 0.5*dt*y2, tk);
        y4 = f(uk + dt*y3, tk);
        
        ukp1 = uk + (1/6)*dt*(y1 + 2*y2 + 2*y3 + y4);
        
        uk = ukp1;
        tk = tkp1;
        
        if min(abs( tkp1-tsv ) ) < 1e-8
            
            u(:,cnt) = uk(1:n-1);
            cnt = cnt + 1;
            
        end
        
    end
    %---
    err(j) = norm( ukp1(1:n-1) - uex(xj(2:end-1),T) ) / norm( uex(xj(2:end-1),T) );
end
toc

%Error plots
figure(3)
c = err(end)*(1./dt^4);
loglog( dtvect, c*(dtvect).^4, 'k--', 'linewidth', 2 ), hold on

%plot err
loglog(dtvect, err , 'b.', 'markersize', 26 )
xlim([1e-2 100])
title('Temporal Convergence RK4')

%make pretty
h = legend('$O(\Delta t^2)$', '$Error$');
set(h, 'Interpreter','latex', 'fontsize', 16, 'Location', 'NorthWest' )
xlabel( '$\Delta t$', 'interpreter', 'latex', 'fontsize', 16)
ylabel( '$||\textbf{e}||/||\textbf{u}_e||$ ', 'interpreter', 'latex', 'fontsize', 16)
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16 )   