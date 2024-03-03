%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical simulation of the Generalized Modified BBM equation
% using RK4
% 
% Periodic boundaries
% 
% Vilma Kjelldahl
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global D1 D3 split_approx

%plotting on/off    1/0
plotting = 0;

%equation, 1-BBM equation. 2-Modified BBM equation
n = 2;

%if n == 2, define which split approximation to use
%1 - eq 14,  2 - eq 16
split_approx = 2;

%number of spatial grid points
mx = 50; 

%order of spatial discretization. 2, 4 or 6
order = 6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%












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


%Time discretization
T = 1.5;     %end time
ht = 0.3 * hx^3; 
mt = floor(T/ht);
tvec = ht:ht:T;
t = 0;


%SBP operators, space
if order == 2
    [D1, D3] = Periodic_C2(L, mx);
elseif order == 4
    [D1, D3] = Periodic_C4(L, mx);
elseif order == 6
    [D1, D3] = Periodic_C6(L, mx);
end
D1 = sparse(D1);
D3 = sparse(D3);



%initial condition, analytic solution
u = BBM_Analytic(0, xvec, L, a,n, c);

% plot
if plotting == 1
    hold on
    p = plot(xvec, u, 'k');
    if n == 1
        ylim([-10 35])
    elseif n == 2
        ylim([-5 10])
    end 
    p.XDataSource = 'xvec';
    p.YDataSource = 'u';
    title(['t = ' num2str(t)])
end


tic %start timing

% RK4 algorithm
for i=1:mt
    w1=RHS(u, a, n);
    w2=RHS(u+ht/2*w1, a, n);
    w3=RHS(u+ht/2*w2, a, n);
    w4=RHS(u+ht*w3, a, n);
    u = u + ht*(w1+2*w2+2*w3+w4)/6;
    t = t+ ht;

    % update plot and print the time every xth iteration. 
    if plotting == 1
        if mod(i, 10) == 0    
            disp(t)
            refreshdata
            drawnow
            title(['t = ' num2str(t)])
        end
    end
end
toc   %end timing


%exact solution
u_exact = BBM_Analytic(t, xvec, L, a, n, c);

%error
error =  l2_norm(u_exact-u, hx);

disp(['l2 error: ' num2str(error)])






%l2 norm
function[e] = l2_norm(vec, h)
    e = sqrt(h)*sqrt(sum(vec.^2));
end

% Analytic solutions - periodic boundaries
function [ u ] = BBM_Analytic(tt, x, L, a, n, c )
    % Analytic solutions to the generalized modified BBM Eq.
    % New travelling wave solutions of different physical structures to generalized BBM equation
    % Physics Letters A 355 (2006) 358–362. Eq. (18)  
    t = rem(tt,L/c);
    um=((n+2)*(n+1)*(c-1)./(a * (1+ cosh(n*sqrt(c-1)*((x+L)-c*t))) )).^(1/n);
    u0=((n+2)*(n+1)*(c-1)./(a * (1+ cosh(n*sqrt(c-1)*((x+0)-c*t))) )).^(1/n);
    up=((n+2)*(n+1)*(c-1)./(a * (1+ cosh(n*sqrt(c-1)*((x-L)-c*t))) )).^(1/n);
    u=um+u0+up;
    u = u';
end

function Vout = RHS(v, a, n)
    global D1 D3 split_approx
    if n == 1
        Vout = -D1*v - a/3*D1*v.^2 - a/3 * v .* D1 *v - D3*v ;
    elseif n == 2
        if split_approx == 1
            Vout = -D1*v - a/4*D1*v.^3 - a/4* v.^2 .* D1*v - D3*v ;
        elseif split_approx == 2
            Vout = -D1*v - a/2* v.*D1*v.^2 - D3*v ;
        end 
    end
end





