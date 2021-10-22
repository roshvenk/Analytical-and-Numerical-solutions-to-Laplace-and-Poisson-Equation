countArray1 = zeros(1,8);
countArray2 = zeros(1,8);

    ndiv = 64;
    nx =  64 + 1;
    nz = 64+ 1;
    L = 25;
    D = 25;
    t0 = 20;
    t1 = 380;
    t2 = 205;


    
    % Create domain
    x = linspace(0, L, nx);
    z = linspace(0, D, nz);
    [X, Z] = meshgrid(x,z);
    
    
    % Initialise solution array for timestep n+1
    % Note indexing: rows = y, columns = x to agree with meshgrid output
    Tnp1 = zeros(nz, nx);
    SOR_Tnp1 = zeros(nz, nx);
    
    
    % Apply boundary conditions. 
    for u = 1:length(x)
        Tnp1(end,u) = t0+ t1*sin(pi*x(u)/L) + t2*sin(5*pi*x(u)/L);      % Bottom boundary.
        SOR_Tnp1(end,u) = t0+ t1*sin(pi*x(u)/L) + t2*sin(5*pi*x(u)/L);  % Bottom boundary.
    end
    Tnp1(:,1) = t0;             % Left boundary.
    Tnp1(:,end) = t0;           % Right boundary.
    Tnp1(1,:) = t0;             % Top boundary.
    SOR_Tnp1(:,1) = t0;         % Left boundary.
    SOR_Tnp1(:,end) = t0;       % Right boundary.
    SOR_Tnp1(1,:) = t0;         % Top boundary.
    
    
    % Initialise error and set tolerance for convergence
    err_GS = 1;
    err_SOR = 1;
    tol = 1e-8;
    count2 = 0;
    
    %% Successive Over-Relaxation Loop
    while  err_SOR > tol
        
        % Update solution array for this timestep
   
        SOR_Tn = SOR_Tnp1;
        
        % Loop over internal points
        for i = 2:nx-1
            for j = 2:nz-1
                
                % SOR
                
                SOR_Tnp1(j,i) = SOR_Tn(j,i) + 1.85*( 0.25*( SOR_Tn(j,i+1) + SOR_Tnp1(j,i-1) + SOR_Tn(j+1,i) + SOR_Tnp1(j-1,i)) - SOR_Tn(j,i));
            end
        end
        
        % Compute error as maximum change in domain (absolute value)
        err_SOR = max(abs(SOR_Tnp1(:) - SOR_Tn(:)));
        count2 = count2 + 1;
        
    end
    countArray2(k) = count2;    

surf(X,Z,SOR_Tnp1,'EdgeColor','none');
view(0, -90);
title('heat map with grid size 64');
xlabel('Length, meters');
ylabel('Depth, meters');
zlabel('Temperature, degrees celsius');