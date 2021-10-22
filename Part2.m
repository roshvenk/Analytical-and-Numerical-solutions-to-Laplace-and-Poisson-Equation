
    nx =  64 + 1;
    nz = 64+ 1;
    L = 25;
    D = 25;
    t0 = 20;
    t1 = 380;
    t2 = 205;


    SORerror = 1;
    tolerance = 1e-8;
    x = linspace(0, L, nx);
    z = linspace(0, D, nz);
    [X, Z] = meshgrid(x,z);
    Tnp1 = zeros(nz, nx);
    Tnp1SOR = zeros(nz, nx);
    
    Tnp1SOR(1,:) = t0;
    Tnp1SOR(:,1) = t0;         
    Tnp1SOR(:,end) = t0;       
    Tnp1(:,1) = t0;             
    Tnp1(:,end) = t0;          
    Tnp1(1,:) = t0;            
    
    for i = 1:length(x)
        Tnp1(end,i) = t0+ t1*sin(pi*x(i)/L) + t2*sin(5*pi*x(i)/L);      
        Tnp1SOR(end,i) = t0+ t1*sin(pi*x(i)/L) + t2*sin(5*pi*x(i)/L);  
    end
            
    while  SORerror > tolerance
        
     
   
        TnSOR = Tnp1SOR;
        
      
        for i = 2:nx-1
            for j = 2:nz-1
                
                
                Tnp1SOR(j,i) = TnSOR(j,i) + 1.85*( 0.25*( TnSOR(j,i+1) + Tnp1SOR(j,i-1) + TnSOR(j+1,i) + Tnp1SOR(j-1,i)) - TnSOR(j,i));
            end
        end
        
        % Compute error as maximum change in domain (absolute value)
        SORerror = max(abs(Tnp1SOR(:) - TnSOR(:)));
    end
      

surf(X,Z,Tnp1SOR,'EdgeColor','none');
view(0, -90);
title('heat map with grid size 64');
xlabel('Length, meters');
ylabel('Depth, meters');
zlabel('Temperature, degrees celsius');