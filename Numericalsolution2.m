%% 2.2: Implementation of the numerical solution using Gauss-Seidel and SOR method. Graphing number of iterations vs each method.
%% Roshan Venkatesan 500316319

countArray1 = zeros(1,8);
countArray2 = zeros(1,8);
nArray = [4,8,16,32,64,128,256,512];
j = 1;
for n = 2:9
    
    %Declaring our constants/points.
    nx = 2^n + 1;
    nz = 2^n + 1;
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
    
    % Apply boundary conditions. 
    for i = 1:length(x)
        Tnp1(end,i) = t0 + t1*sin(pi*x(i)/L) + t2*sin(5*pi*x(i)/L);      % Bottom boundary.
        Tnp1SOR(end,i) = t0 + t1*sin(pi*x(i)/L) + t2*sin(5*pi*x(i)/L);  % Bottom boundary.
    end
    Tnp1(:,1) = t0;             % Left boundary.
    Tnp1(:,end) = t0;           % Right boundary.
    Tnp1(1,:) = t0;             % Top boundary.
    Tnp1SOR(:,1) = t0;         % Left boundary.
    Tnp1SOR(:,end) = t0;       % Right boundary.
    Tnp1SOR(1,:) = t0;         % Top boundary.
    
    % Initialise error and set tolerance for convergence
    GSerror = 1;
    SORerror = 1;
    t = 1e-8;      %Tolerance value.

    [count2] = SOR(SORerror,t,nx,nz,Tnp1SOR);
    %Add our iteration count to this array for SOR.
    countArray2(j) = count2;

    
    [count1] = GS(GSerror,t,nx,nz,Tnp1);
    
    %add our iteration count to this array for G-S.
    countArray1(j) = count1;
    j = j + 1;
end

%We remove all the zeroes from our iteration arrays
countArray1 = countArray1(countArray1~=0);
countArray2 = countArray2(countArray2~=0);

%Plot our countArray vs the grid sizes.
hold on;
grid on;
ylim([0 250000]);
xlim([0 600]);
plot(nArray,countArray1);
plot(nArray,countArray2);
ylabel('Number of iterations');
xlabel('Grid size');
legend('Gauss method', 'SOR Method');

function [count] = SOR(SORerror, tolerance, nx,nz, Tnp1SOR)
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

function [count] = GS(GSerror, tolerance, nx,nz, Tnp1)
count = 0;
%Our G-S loop
    while  GSerror > tolerance
        
        %In each timestep, we update solution.
        Tn = Tnp1;
        
        for i = 2:nx-1
            for j = 2:nz-1
                
                % Gauss Seidel
                Tnp1(j,i) = (1/4) * (Tn(j,i+1) + Tnp1(j,i-1) + Tn(j+1,i) + Tnp1(j-1,i));   
                
            end
        end
    
        % Update our error for Gauss-Seidel.
        GSerror = max(abs(Tnp1(:) - Tn(:)));
        count = count + 1;      
        
    end
end