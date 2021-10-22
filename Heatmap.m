%% 2.3 Heatmap for grid space of 2^8
countArray2 = zeros(1,8);


    %Declaring our constants/points.
    nx = 2^3 + 1;
    nz = 2^3 + 1;
    L = 25;
    D = 25;
    t0 = 20;
    t1 = 380;
    t2 = 205;

    % Grid Spacing.
    dx = L/(nx-1);
    dz = D/(nz-1);
    
    % Create domains in the x and z direction.
    x = linspace(0, L, nx);
    z = linspace(0, D, nz);
    [X, Z] = meshgrid(x,z);
    
    
    %Initialize our solution arrays.
    Tnp1 = zeros(nz, nx);
    Tnp1SOR = zeros(nz, nx);
    
    % Initialize boundary conditions. 
    for i = 1:length(x)
        Tnp1(end,i) = t0 + t1*sin(pi*x(i)/L) + t2*sin(5*pi*x(i)/L);     
        Tnp1SOR(end,i) = t0 + t1*sin(pi*x(i)/L) + t2*sin(5*pi*x(i)/L);  
    end
    Tnp1(:,1) = t0;       
    Tnp1(:,end) = t0;           
    Tnp1(1,:) = t0;             
    Tnp1SOR(:,1) = t0;         
    Tnp1SOR(:,end) = t0;       
    Tnp1SOR(1,:) = t0;        
    
    % Initialise error and set tolerance for convergence
    SORerror = 1;
    t = 1e-8;      %Tolerance value.
    
    
     
    %Add our iteration count to this array for SOR.
    Tnp1SOR = SOR(SORerror, t, nx,nz,Tnp1SOR);

    
%Plot our heatmap for grid-size 8
surf(X,Z,Tnp1SOR,'EdgeColor','none');
view(0, -90);
title('heat map with size 2^8');
xlabel('Length, meters');
ylabel('Depth, meters');
zlabel('Temperature, degrees celsius');

function Tnp1SOR = SOR(SORerror, tolerance, nx,nz, Tnp1SOR)
count = 0;
while  SORerror > tolerance
        
        % Update solution array for this timestep
        
        TnSOR = Tnp1SOR;
        
        
        for i = 2:nx-1
            for j = 2:nz-1
                
                % SOR
                Tnp1SOR(j,i) = TnSOR(j,i) + 1.85*( 0.25*( TnSOR(j,i+1) + Tnp1SOR(j,i-1) + TnSOR(j+1,i) + Tnp1SOR(j-1,i) ) - TnSOR(j,i) );
                
            end
        end
        
        %Update our error for SOR.
        SORerror = max(abs(Tnp1SOR(:) - TnSOR(:)));
        %Increment our iterations.
        count = count + 1;
        
end
end
