%% 2.4 Error norms for SOR Method.
nArray = [4,8,16,32,64,128,256,512];

%Create our empty arrays for our error norms calculations.
error1 = zeros(1,length(nArray));
error2 = zeros(1,length(nArray));
errorInf = zeros(1,length(nArray));
order1 = zeros(1,length(nArray));
order2 = zeros(1,length(nArray));
orderinf = zeros(1,length(nArray));

k = 1;
for n = 2:9
    nx = 2^n + 1;
    nz = 2^n + 1;
    L = 25;
    D = 25;
    t0 = 20;
    t1 = 380;
    t2 = 205;

    % Create domain
    x = linspace(0, L, nx);
    z = linspace(0, D, nz);
    [X, Z] = meshgrid(x,z);

    %Solution domains for numerical and analytical
    Tnp1SOR = zeros(nz, nx);
    Tanalytical = zeros(nz,nx);
    
    % boundary conditions. 
    for i = 1:length(x)      
        Tnp1SOR(end,i) = t0+ t1*sin(pi*x(i)/L) + t2*sin(5*pi*x(i)/L); 
    end           
    Tnp1SOR(:,1) = t0;        
    Tnp1SOR(:,end) = t0;       
    Tnp1SOR(1,:) = t0;         
    
     % Grid Spacing.
    dx = L/(nx-1);
    dz = D/(nz-1);
    
    % Initialise error and set tolerance for convergence
    SORerror = 1;
    t = 1e-8;    
    Tnp1SOR = SOR(SORerror,t,nx,nz,Tnp1SOR);
    Tanalytical = analyticalSolution(Tanalytical,nx,nz,t0,t1,t2,x,z,D,L);
          
            % Go through each point in the domain
         
       
         
         %Calculating the difference between our analytical and numerical
         %Reshaping it to prepare for norms calculation
            diff = Tnp1SOR(2:end-1, 2:end-1) - Tanalytical(2:end-1, 2:end-1);
            diffVec = reshape(diff',1,numel(diff));
            error1(k) = norm(diffVec,1)*dx*dz;
            error2(k) = norm(diffVec,2)*sqrt(dx*dz);
            errorInf(k) = norm(diffVec,Inf);
            
            if k > 1
            order1(k) = (log(error1(k-1)/error1(k)))/(log(nArray(k)/nArray(k-1)));
            order2(k) = (log(error2(k-1)/error2(k)))/(log(nArray(k)/nArray(k-1)));
            orderinf(k) = (log(errorInf(k-1)/errorInf(k)))/(log(nArray(k)/nArray(k-1)));
            end
        
        k = k + 1;
end

T = table(nArray',error1',error2',errorInf',order1',order2',orderinf','VariableNames',["N","L1","L2","LInf","O1","O2","OInf"]);


%% Function to calculate analytical solution.
function Tanalytical = analyticalSolution(Tanalytical,nx,nz,t0,t1,t2,x,z,D,L)

        for i = 1:nx
             
          for j = 1:nz
             
            % Sum all terms in the homogeneous solution at point
            Tanalytical(j,i) = t0 + (t1/sinh(pi*D/L))*sin(pi*x(i)/L)*sinh(pi*z(j)/L) + (t2/sinh(5*pi*D/L)) *sin(5*pi*x(i)/L)*sinh(5*pi*z(j)/L);
          end
          
        end
end

%% Function to calculate using SOR method.
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
