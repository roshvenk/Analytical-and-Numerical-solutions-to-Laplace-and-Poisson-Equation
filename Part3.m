% Gas constant (J kg^-1 K^-1)
r = 287.058; 
% Gravitational acceleration (ms^-2) 
g = 9.81;
    L=25;
    D=25;
    t0 = 20;
    t1 = 380;
    t2 = 205;
    
    %Grid size
    ndiv = 256;
    nx = ndiv + 1;
    nz = ndiv + 1;
    
    % Create domain
    x = linspace(0, L, nx);
    z = linspace(0, D, nz);
    % Calculate grid spacing and maximum stable timestep
    dx = L/(nx-1);
    dz = D/(nz-1);
    
    
    % Our temperature field
    Tanalytical = zeros(nz, nx);           
    for i = 1:nx
        for j = 1:nz        
            Tanalytical(j,i) = t0 + (t1/sinh(pi*D/L))*sin(pi*x(i)/L)*sinh(pi*z(j)/L) + (t2/sinh(5*pi*D/L)) *sin(5*pi*x(i)/L)*sinh(5*pi*z(j)/L);  
        end
    end
    % Convert to Kelvin
    Tanalytical = Tanalytical + 273.15; 
    
    % Initialise Pressure arrays where every point is 101325. 
    Pressure = ones(nz,nx) * 101325;       
    Pressurep1 = ones(nz,nx) * 101325;      
    
    
    %Density array calculation
    Density = Pressure./(Tanalytical * r);                                     
    InitialD = Density;
    
    
    % Initialise error and set tolerance for convergence
    S3_err_Pn = 1;
    S3_tol = 1e-8;
    
    while S3_err_Pn > S3_tol
        
        % Update our pressure solution for this iteration
        Pressure = Pressurep1;        
        
        
        % Calculating source term and find new pressure field on it.
        for i = 2:nx-1
            for j = 2:nz-1
                pressuredz = (Pressure(j+1,i)-Pressure(j-1,i))/(2*dz);
                densitydz = (Density(j+1,i)-Density(j-1,i))/(2*dz);
                pressuredx = (Pressure(j,i+1)-Pressure(j,i-1))/(2*dx);
                densitydx = (Density(j,i+1)-Density(j,i-1))/(2*dx);
                source = ((-dz^2)/Density(j,i))*((pressuredz*densitydz) + (pressuredx*densitydx));
                
                % Pressure solution.
                Pressurep1(j,i) = Pressure(j,i) + 1.85*( source + 0.25*( Pressure(j,i+1) + Pressurep1(j,i-1) + Pressure(j+1,i) + Pressurep1(j-1,i) ) - Pressure(j,i) );
            end
        end
       
       % Calculating new density.
        for i = 1:nx
            for j = 1:nz
                
                % Ideal Gas Law
                Density(j,i) = Pressurep1(j,i)/(Tanalytical(j,i) * r);
                
            end
        end
        
            % Top boundary
            Pressurep1(1,:) = 101325;
            % Left boundary
            Pressurep1(:,1) = Pressurep1(:,2);
            % Right boundary
            Pressurep1(:,nx) = Pressurep1(:,nx-1);
            
            
          for i = 1:nz
            Pressurep1(end,i) = Pressurep1(end-1,i) + Density(end-1,i)*g*dz;
          end
       
        
        
        % Part (vi) -----------------------------------------------------
        Density(:,1) = Pressurep1(:,1)./(Tanalytical(:,1) * r);
        Density(1,:) = Pressurep1(1,:)./(Tanalytical(1,:) * r);
        Density(end,:) = Pressurep1(end,:)./(Tanalytical(end,:) * r);
        Density(:,end) = Pressurep1(:,end)./(Tanalytical(:,end) * r);

      % Compute error as maximum change in domain (absolute value).
        S3_err_Pn = max(abs(Pressurep1(:) - Pressure(:)));
        
    end
      
figure()
surf(x,z,Density),shading interp,title('Density'),xlabel('x'),ylabel('y'),colorbar;
view(0,-90);
figure()
surf(x,z,Pressurep1),shading interp,title('Pressure'),xlabel('x'),ylabel('y'),colorbar;
view(0,-90);
figure()
surf(x,z,InitialD),shading interp,title('Initial Density'),xlabel('x'),ylabel('y'),colorbar;
view(0,-90);