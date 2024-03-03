function [D1, HI] = SBP4(m,h)
% Create diagonal norm FD-SBP operators for the first derivative
% Fourth order accurate in the interior and second order accurate at the
% boundaries
% 8 boundary points

        H=diag(ones(m,1),0);
        H(1:4,1:4)=diag([17/48 59/48 43/48 49/48]);
        H(m-3:m,m-3:m)=fliplr(flipud(diag([17/48 59/48 43/48 49/48])));
        HI = inv(H*h);
          
        D1=(-1/12*diag(ones(m-2,1),2)+8/12*diag(ones(m-1,1),1)- ...
            8/12*diag(ones(m-1,1),-1)+1/12*diag(ones(m-2,1),-2));
        
        D1(1:4,1:6)=[-24/17,59/34,-4/17,-3/34,0,0; -1/2,0,1/2,0,0,0; 4/43,-59/86,0,59/86,-4/43,0; 3/98,0,-59/98,0,32/49,-4/49];
        D1(m-3:m,m-5:m)=flipud( fliplr(-D1(1:4,1:6)));
        D1=D1/h;
    end