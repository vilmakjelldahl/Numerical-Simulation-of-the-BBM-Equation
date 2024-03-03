function [D1, D3, D5 ] = Periodic_C6( L,n )
% Create periodic FD opereators for first,
% third and fifth derivative. Sixth order accurate.
% On a periodic domain of width L=(b-a), using n points.
% The x points are given by x(i) = a+i*h, i=0,1,...n-1
% where h=L/n and f(a)=f(b) due to periodicity. 

%Author: Ken Mattsson


h=L/n;

% D1 operator
l=3;
r=3;
d1=[-1/60 3/20 -3/4 0 3/4 -3/20 1/60];
v=zeros(1,n);
for i=1:r+1
    v(i)=d1(i+l);
end
for i=1:l
    v(n-i+1)=d1(l-i+1);
end

D1=1/h*toeplitz(circshift(flipud(v(:)),1),v);


% D3 operator
l=4;
r=4;
d3=[-7/240 3/10 -169/120 61/30 0 -61/30 169/120 -3/10 7/240];
v=zeros(1,n);
for i=1:r+1
    v(i)=d3(i+l);
end
for i=1:l
    v(n-i+1)=d3(l-i+1);
end

D3=1/h^3*toeplitz(circshift(flipud(v(:)),1),v);

% D5 operator
l=5;
r=5;
d5=[-13/288 19/36 -87/32 13/2 -323/48 0 323/48 -13/2 87/32 -19/36 13/288];
v=zeros(1,n);
for i=1:r+1
    v(i)=d5(i+l);
end
for i=1:l
    v(n-i+1)=d5(l-i+1);
end

D5=1/h^5*toeplitz(circshift(flipud(v(:)),1),v);



end

