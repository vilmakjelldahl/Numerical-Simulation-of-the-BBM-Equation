function [D1, D3, D5 ] = Periodic_C2( L,n )
% Create periodic FD opereators for first,
% third and fifth derivative. Second order accurate.
% On a periodic domain of width L=(b-a), using n points.
% The x points are given by x(i) = a+i*h, i=0,1,...n-1
% where h=L/n and f(a)=f(b) due to periodicity. 

%Author: Ken Mattsson

h=L/n;

% D1 operator
l=1;
r=1;
d1=[-1/2 0 1/2];
v=zeros(1,n);
for i=1:r+1
    v(i)=d1(i+l);
end
for i=1:l
    v(n-i+1)=d1(l-i+1);
end
D1=1/h*toeplitz(circshift(flipud(v(:)),1),v);


% D3 operator
l=2;
r=2;
d3=[-1/2 1 0 -1 1/2];
v=zeros(1,n);
for i=1:r+1
    v(i)=d3(i+l);
end
for i=1:l
    v(n-i+1)=d3(l-i+1);
end
D3=1/h^3*toeplitz(circshift(flipud(v(:)),1),v);

% D5 operator
l=3;
r=3;
d5=[-1/2 2 -5/2 0 5/2 -2 1/2];
v=zeros(1,n);
for i=1:r+1
    v(i)=d5(i+l);
end
for i=1:l
    v(n-i+1)=d5(l-i+1);
end
D5=1/h^5*toeplitz(circshift(flipud(v(:)),1),v);

end

