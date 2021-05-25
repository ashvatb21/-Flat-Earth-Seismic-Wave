%% Trapezoidal

clear all, close all, clc

% parameters for problem (Distance Reduced by a Factor of 10)
a = 0; % Lower Bound x (km)
b = 40075/10; % Upper Bound x (km)
c = 5.86/10; % Wave Speed (km/s)
A = 0.002/10; % Amplitude (km)
ln = b-a; % The domian of x (km)
T = 60*60*3; % Final time to run to (seconds)
dt = 1; % Make dt small in spatial convergence test so that time error doesn't pollute convergence


% Initial Conditions, Gaussian Function
std = 87.5; % Roundness of the waveform
u_init = @(x) A*exp(-((x-(b-(9250/10)))/(2*std)).^2); % Initial Displacement
v_init = @(x) 0.*x; % Initial Velocity

%g(x,t) and Boundary Conditions
bc_a = @(t) 0; % Dirichlet Boundary Condition at 'a'
bc_b = @(t) 0; % Dirichlet Boundary Condition at 'b'
fcn = @(x,t) zeros(length(x),1); % Forcing Function

%# of n points to use
nvect = [200];

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
    
end
toc

[X,T] = meshgrid( xj(2:end-1), tsv );
U = u(1:n-1,:).';

%--Waterfall plot of solns
figure(1)
waterfall( X,T,U ), hold on
set( gca, 'fontsize', 15, 'ticklabelinterpreter', 'latex' )
title('$u_{FD}$', 'fontsize', 20, 'interpreter', 'latex')
xlabel('$x$', 'fontsize', 15, 'interpreter', 'latex')
ylabel('$t$', 'fontsize', 15, 'interpreter', 'latex')
hold off

% Movie for the wave travel
figure(2)
for j = 1:nsnps
    plot(X(1,:),U(j,:));
    xlim([a b]);ylim([-0.04/10 0.04/10]);
    F(j) = getframe(gcf) ;
    drawnow
end

% save video
writerObj = VideoWriter('SeismicWave.avi');
% set the seconds per image
writerObj.FrameRate = 10;
% open the video writer
open(writerObj);

% write the frames to the video
for j = 1:length(F)
    % convert the image to a frame
    frame = F(j) ;
    writeVideo(writerObj, frame);
end

% close the writer object
close(writerObj);

% Snapshots of the wave at different times
figure(3)
subplot(4,1,1); hold on;
H = area(X(1,:),U(1,:));
set(H,'FaceColor',[1 0 0]);
xlim([a b]);
ylim([-0.04/100 0.04/100]);
xlabel('Distance in Kilometers');
title('Wave Snapshots at different time intervals');

subplot(4,1,2);
H = area(X(1,:),U(83,:));
set(H,'FaceColor',[1 0 0]);
xlim([a b]);
ylim([-0.04/100 0.04/100]);
xlabel('Distance in Kilometers');
% title('Wave Snapshots at different intervals');

subplot(4,1,3);
H = area(X(1,:),U(167,:));
set(H,'FaceColor',[1 0 0]);
xlim([a b]);
ylim([-0.04/100 0.04/100]);
xlabel('Distance in Kilometers');
% title('Wave Snapshots at different intervals');

subplot(4,1,4);
H = area(X(1,:),U(250,:));
set(H,'FaceColor',[1 0 0]);
xlim([a b]);
ylim([-0.04/100 0.04/100]);
xlabel('Distance in Kilometers');
% title('Wave Snapshots at different intervals');

hold off

%% RK4

clear all; close all; clc

% parameters for problem (Distance Reduced by a Factor of 10)
a = 0; % Lower Bound x (km) (Distance Reduced by a Factor of 10)
b = 40075/10; % Upper Bound x (km)
c = 5.86/10; % Wave Speed (km/s)
A = 0.002/10; % Amplitude (km)
ln = b-a; % The domian of x (km)
T = 60*60*3; % Final time to run to (seconds)
dt = 1; % Make dt small in spatial convergence test so that time error doesn't pollute convergence


% Initial Conditions, Gaussian Function
std = 87.5; % Roundness of the waveform
% u_init = @(x) A*exp(-((x-(a+b)/2)/(std)).^2); % Initial Displacement
u_init = @(x) A*exp(-((x-(b-(9250/10)))/(2*std)).^2); % Initial Displacement
v_init = @(x) 0.*x; % Initial Velocity

%g(x,t) and Boundary Conditions
bc_a = @(t) 0; % Dirichlet Boundary Condition at 'a'
bc_b = @(t) 0; % Dirichlet Boundary Condition at 'b'
fcn = @(x,t) zeros(length(x),1); % Forcing Function

%# of n points to use
nvect = [200];

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
        
        %It sets things up to only save for a relatively
        %small # of times.
        if min(abs( tkp1-tsv ) ) < 1e-8
            
            u(:,cnt) = uk(1:n-1);
            cnt = cnt + 1;
            
        end
        
    end
end
toc

[X,T] = meshgrid( xj(2:end-1), tsv );

%--Waterfall plot of solns
figure(1)
waterfall( X,T, u.' ), hold on
set( gca, 'fontsize', 15, 'ticklabelinterpreter', 'latex' )
title('$u_{FD}$', 'fontsize', 20, 'interpreter', 'latex')
xlabel('$x$', 'fontsize', 15, 'interpreter', 'latex')
ylabel('$t$', 'fontsize', 15, 'interpreter', 'latex')
hold off

% Movie for the wave travel
figure(2)
for j = 1:nsnps
    U = u.';
    plot(X(1,:),U(j,:));
    xlim([a b]);ylim([-0.04/10 0.04/10]);
    drawnow
end

% Snapshots of the wave at different times
figure(3)
subplot(4,1,1); hold on;
H = area(X(1,:),U(1,:));
set(H,'FaceColor',[1 0 0]);
xlim([a b]);
ylim([-0.04/100 0.04/100]);
xlabel('Distance in Kilometers');
title('Wave Snapshots at different time intervals');

subplot(4,1,2);
H = area(X(1,:),U(83,:));
set(H,'FaceColor',[1 0 0]);
xlim([a b]);
ylim([-0.04/100 0.04/100]);
xlabel('Distance in Kilometers');
% title('Wave Snapshots at different intervals');

subplot(4,1,3);
H = area(X(1,:),U(167,:));
set(H,'FaceColor',[1 0 0]);
xlim([a b]);
ylim([-0.04/100 0.04/100]);
xlabel('Distance in Kilometers');
% title('Wave Snapshots at different intervals');

subplot(4,1,4);
H = area(X(1,:),U(250,:));
set(H,'FaceColor',[1 0 0]);
xlim([a b]);
ylim([-0.04/100 0.04/100]);
xlabel('Distance in Kilometers');
% title('Wave Snapshots at different intervals');