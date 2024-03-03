%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical simulation of the Generalized Modified BBM equation
% using SBP-SAT in time
% 
% Periodic boundaries
% 
% Vilma Kjelldahl
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%global variables
global a n c D1x It A R u0 D1t Ix mx mt split_approx

%plot on/off  1/0
plotting = 0;

%quadrature, 1-Gauss-Lobatto. 2-SBP(4,2) 3-SBP(6,3)
quadrature = 1;

%newton tolerance: set to the expected order of accuracy or lower
tol = 1e-8;

%equation, 1-BBM equation. 2-Modified BBM equation
n = 2;

%if n == 2, define which split approximation to use
%1 - eq 28,  2 - eq 41
split_approx = 2;

%number of spatial grid points
mx = 100;     

%order of spatial discretization. 2, 4 or 6
order = 6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5







%equation parameters
a = 1;
c = 10;

% Space discretization
%Periodic BC, rightmost grid point left out
xl = -5;    %left boundary
xr = 5;     %right boundary
L = xr-xl;
hx = L/mx;  
xvec = linspace(xl, xr-hx, mx);  


%SBP operators space, periodic
if order == 2
    [D1x, D3x] = Periodic_C2(L, mx);
elseif order == 4
    [D1x, D3x] = Periodic_C4(L, mx);
elseif order == 6
    [D1x, D3x] = Periodic_C6(L, mx);
end
D1x = sparse(D1x);
D3x = sparse(D3x);
Ix = sparse(eye(mx));


%Time discretization
sigma = -1;
T = 1.5;          %end time
%temporal resolution
if quadrature == 1
    ht = 0.1*hx;     
elseif quadrature == 2
    ht = 0.04*hx;    
elseif quadrature == 3
    ht = 0.1*hx;    
end

nt = round(T/ht);
n_blocks = nt;
block = T/n_blocks;    %time block size



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%quadrature specific operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if quadrature == 1
    %Gauss lobatto
    mt = 4;
    It = eye(mt);
    
    T_GL = [-1 -1/5*sqrt(5) 1/5*sqrt(5) 1];
    
    e_1t=zeros(mt,1);e_1t(1) = 1;
    e_mt=zeros(mt,1);e_mt(mt) = 1;
    E0 = diag(e_1t);
    
    Ht = eye(mt)*5; 
    Ht(1,1) = 1; Ht(end,end) = 1; 
    Ht = (0.5*block)/6 * Ht;
    HIt = inv(Ht);
    
    D1t  =1/(0.5*block)*...
             [-3                   -(5*sqrt(5))/(sqrt(5)-5) -(5*sqrt(5))/(sqrt(5)+5)   1/2;...
            (sqrt(5))/(sqrt(5)-5)   0                      sqrt(5)/2               -(sqrt(5))/(sqrt(5)+5);...
            (sqrt(5))/(sqrt(5)+5)  -sqrt(5)/2               0                   -(sqrt(5))/(sqrt(5)-5);...
            -1/2                 (5*sqrt(5))/(sqrt(5)+5) (5*sqrt(5))/(sqrt(5)-5)            3];
    
elseif quadrature == 2
    %SBP4
    mt = 8;
    It = eye(mt);
    ht_i = block/(mt-1);
    
    [D1t, HIt] = SBP4(mt, ht_i);

    e_1t=zeros(mt,1);e_1t(1) = 1;
    e_mt=zeros(mt,1);e_mt(mt) = 1;
    E0 = diag(e_1t);

elseif quadrature == 3
    %SBP6
    mt = 12;
    It = eye(mt);
    ht_i = block/(mt-1);
    
    [D1t, HIt] = SBP6(mt, ht_i);

    e_1t=zeros(mt,1);e_1t(1) = 1;
    e_mt=zeros(mt,1);e_mt(mt) = 1;
    E0 = diag(e_1t);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dt_bar = sparse(kron((D1t - sigma*HIt*E0), Ix)); 
Dx_bar = sparse(kron(It, -D1x - D3x));
R = sparse( kron(diag(sigma*HIt*e_1t), Ix) );

A = Dt_bar - Dx_bar;
[L_, U_, P] = lu(A);   %LU-factorizing






%Initial condition
f = BBM_Analytic(0, xvec, L);
ONE = ones(mt,1);
e_end = kron(e_mt', Ix);

U = [f];   %matrix saving the solutions in each time-step, used for plotting
tic
for i_block = 1:n_blocks
    du = 1;
    u = kron(ONE,f);
    u0 = u;

    while abs(normest(du)) > tol
        %construct b
        b = (Dx_bar-Dt_bar)*u - R*u0 + F(u);

        %solve A*du = b
        y = L_\(P*b);
        du = U_\y;

        %update
        u = u+du;
    end
    f = e_end*u;

    U = horzcat(U,f);   %save solution 
end
toc


%compute error
u_exact = BBM_Analytic(T, xvec, L);
error = l2_norm(f-u_exact, hx);

disp(['l2 error: ' num2str(error)])




if plotting == 1
    uu = U(1:mx,1);

    hold on
    if n == 1
        ylim([-10 37])
    elseif n == 2
        ylim([-8 10])
    end
    p = plot(xvec, uu, 'b');
    p.XDataSource = 'xvec';
    p.YDataSource = 'uu';
    
    for i=1:n_blocks
        if mod(i,1) == 0
            uu = U(1:end,i+1);
        
            refreshdata
            drawnow
        end
    end
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[fout] = F(u)
    global a D1x It n split_approx
    D = sparse(kron(It, D1x));
    if n == 1
        fout = -a/3 * (D * u.^2 + u.*D*u);
    elseif n == 2
        if split_approx == 1
            fout = -a/4 *(D*u.^3 + u.^2 .* D*u);
        elseif split_approx == 2
            fout = -a/2* u.*D*u.^2;
        end
    end 
end

function[e] = l2_norm(vec, h)
    e = sqrt(h)*sqrt(sum(vec.^2));
end

%analytic solution, periodic
function [ u ] = BBM_Analytic(tt, x, L)
    % Analytic solutions to the modified BBM Eq.
    % New travelling wave solutions of different physical structures to generalized BBM equation
    % Physics Letters A 355 (2006) 358–362. Eq. (18)
    global a n c
    t = rem(tt,L/c);
    um=((n+2)*(n+1)*(c-1)./(a * (1+ cosh(n*sqrt(c-1)*((x+L)-c*t))) )).^(1/n);
    u0=((n+2)*(n+1)*(c-1)./(a * (1+ cosh(n*sqrt(c-1)*((x+0)-c*t))) )).^(1/n);
    up=((n+2)*(n+1)*(c-1)./(a * (1+ cosh(n*sqrt(c-1)*((x-L)-c*t))) )).^(1/n);
    u=um+u0+up;
    u = u';
end


