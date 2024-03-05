%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical simulation of the Generalized Modified BBM equation
% using RK4
% 
% Non-periodic boundaries imposed with projection
% Dirichlet and neumann on the right boundary. Dirichlet on the left
% boundary.
% 
% Vilma Kjelldahl
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global D1 D3 P HI e_l e_r d1_l d1_r d2_r d2_l split_approx

%plotting on/off    1/0
plotting = 1;

%equation, 1-BBM equation. 2-Modified BBM equation
n = 2;

%if n == 2, define which split approximation to use
%1 - eq 14,  2 - eq 16
split_approx = 1;

%load reference solution
if n == 1
    load('REF_0.1h_proj_-10_n1.mat')
elseif n == 2
    if split_approx == 1
        load('REF_0.1h_proj_-10_n2.mat')
    elseif split_approx == 2
        load('REF_0.1h_proj_-10_n22.mat')
    end
end 

%number of spatial grid points
mx = 51; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%equation parameters
a = 1;
c = 10;

% Space discretization
xl = -5;    %left boundary
xr = 5;     %right boundary
L = xr-xl;
hx = L/(mx-1);  
xvec = linspace(xl, xr, mx);


%Time discretization
T = 0.6;     %end time
ht = 0.3 * hx^3; 
mt = floor(T/ht);
tvec = 0:ht:T;
t = 0;


%SBP operators
order = 6;
[D1, D2, D3, HI, e_l, e_r, d1_l, d1_r, d2_r, d2_l] = SBP_operators(order, mx, hx);
HI = sparse(HI);
D1 = sparse(D1);
D3 = sparse(D3);

Ix = eye(mx);

% BC projection
if n == 1
    L = [e_l' ; e_r' ; d1_r];
    P = Ix - HI*L'*inv(L*HI*L')*L;
elseif n == 2
    L = [e_l' ; e_r' ; d1_r];
    P = Ix - HI*L'*inv(L*HI*L')*L;
end





%initial condition
u = BBM_Analytic(0, xvec, L, a,n, c);

% plot
if plotting == 1
    hold on
    p = plot(xvec, u, 'k');
    if n == 1
        ylim([-20 35])
    elseif n == 2
        ylim([-10 13])
    end 
    p.XDataSource = 'xvec';
    p.YDataSource = 'u';
    title(['t = ' num2str(t)])
end


% RK4 algorithm
tic  %start timing
for i=1:mt
    w1=RHS(u, a, n, c);
    w2=RHS((u+ht/2*w1),  a, n, c);
    w3=RHS( (u+ht/2*w2),  a, n, c);
    w4=RHS((u+ht*w3), a, n, c);
    u = (u + ht*(w1+2*w2+2*w3+w4)/6);
    t = t+ ht;

    % update plot every xth iteration. 
    if plotting == 1
        if mod(i, 10) == 0    
            refreshdata
            drawnow
            title(['n = ' num2str(n) ', t = ' num2str(t)])
        end
    end
end
toc   %end timing



%compute error
error = log10(l2_norm((u-u_exact_approx(1:2000/(mx-1):end)), hx));
disp(['Estimated l2 error: ' num2str(error)])




function[e] = l2_norm(vec, h)
    e = sqrt(h)*sqrt(sum(vec.^2));
end

%Analytic solution, non-periodic
function [ u ] = BBM_Analytic(t, x, L, a, n, c )
    % Analytic solutions to the modified BBM Eq.
    % New travelling wave solutions of different physical structures to generalized BBM equation
    % Physics Letters A 355 (2006) 358–362. Eq. (18)

    u=((n+2)*(n+1)*(c-1)./(a * (1+ cosh(n*sqrt(c-1)*(x-c*t))) )).^(1/n);
    u = u';
end


function Vout = RHS(v, a, n, c)
    global D1 D3 P split_approx
    w = P*v;
    if n == 1
        Vout = -P*D1*w - a/3*P*D1*w.^2 - a/3*w .* P*D1 *w - P*D3*w;
    elseif n == 2
        if split_approx == 1
            Vout = -P*D1*w - a/4*P*D1*w.^3 - a/4*w.^2 .*P* D1*w - P*D3*w ;
        elseif split_approx == 2
            Vout = -P*D1*w - a/2*w.*P*D1*w.^2 - P*D3*w ;
        end
    end
end



