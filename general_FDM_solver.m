function [phi, x, y]=general_FDM_solver(Lx, Ly, Nx, Ny, f1, f2, h1, h2, phi0)

    x=linspace(0, Lx, Nx);
    y=linspace(0, Ly, Ny);
    dx=x(2)-x(1);
    dy=y(2)-y(1);
    A=spalloc(Nx*Ny, Nx*Ny,5*Nx*Ny);
    rhs_val=zeros(Nx*Ny, 1);

    for i=1:Nx % bottom boundary
     k=i;
     A(k,k)=1;
     rhs_val(k,1)=phi0;
    end
    for j=2:Ny % left boundary
     k=(j-1)*Nx+1;
     A(k,k)=1;
     rhs_val(k,1)=phi0;
    end
    for i=2:Nx  % top boundary
        k=(Ny-1)*Nx+i;
        A(k,k)=1+h2*dy;
        A(k,k-Nx)=-1;
    end
    for j=2:Ny-1 %right boundary
        k=j*Nx;
        A(k,k)=1+h1*dx;
        A(k,k-1)=-1;
    end
    for i=2:Nx-1
        for j=2:Ny-1
            k=(j-1)*Nx+i;
            A(k,k)=-2*dy^2/dx^2-2-f1(x(i),y(j))*dy^2;
            A(k,k-1)=dy^2/dx^2;
            A(k,k+1)=dy^2/dx^2;
            A(k,k-Nx)=1;
            A(k,k+Nx)=1;
            rhs_val(k,1)=f2(x(i),y(j))*dy^2;
        end
    end
    phi_vec = A\rhs_val;             
    phi = reshape(phi_vec, Nx, Ny);  