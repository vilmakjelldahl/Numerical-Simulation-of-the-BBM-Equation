%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2nd, 4th and 6th order FD-SBP operators for the 1st, 2nd and 3rd
% derivative. Constructed by Ken Mattsson.
%
% 6 boundary points
% Diagonal norm
% All operators are built on same H-norm, which is necessary if used together.
% 
% 
% 
% H      Norm
% D1     approx 1st derivative
% D2     approx 2nd derivative
% D3     approx 3rd derivative
%
% (m) - number of grid points,  (h) - step length 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function[D1, D2, D3, HI, e_l, e_r, d1_l, d1_r, d2_r, d2_l] = SBP_operators(order, m, h)

    if order == 2
        H=diag(ones(m,1),0);H(1,1)=1/2;H(m,m)=1/2;
        H=H*h;
        HI=inv(H);

        % First derivative SBP operator, 1st order accurate at first 6 boundary points
        q1=1/2;
        Q=q1*(diag(ones(m-1,1),1)-diag(ones(m-1,1),-1));
        e_l=zeros(m,1);e_l(1)=1;
        e_r=zeros(m,1);e_r(m)=1;
        D1=HI*(Q-1/2*e_l*e_l'+1/2*e_r*e_r') ;

        % Second derivative, 1st order accurate at first 6 boundary points
        m1=-1;m0=2;
        M=m1*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1))+m0*diag(ones(m,1),0);M(1,1)=1;M(m,m)=1;
        M=M/h;
        d1_U=[-3/2 2 -1/2]/h;
        d1_l=zeros(1,m);
        d1_l(1:3)=d1_U;
        d1_r=zeros(1,m);
        d1_r(m-2:m)=fliplr(-d1_U);
        D2=HI*(-M-e_l*d1_l+e_r*d1_r);

        % Third derivative, 1/h order accurate at first 6 boundary points
        q2=1/2;q1=-1;
        Q3=q2*(diag(ones(m-2,1),2)-diag(ones(m-2,1),-2))+q1*(diag(ones(m-1,1),1)-diag(ones(m-1,1),-1));   
        Q3_U = [0 -0.13e2 / 0.16e2 0.7e1 / 0.8e1 -0.1e1 / 0.16e2; 0.13e2 / 0.16e2 0 -0.23e2 / 0.16e2 0.5e1 / 0.8e1; -0.7e1 / 0.8e1 0.23e2 / 0.16e2 0 -0.17e2 / 0.16e2; 0.1e1 / 0.16e2 -0.5e1 / 0.8e1 0.17e2 / 0.16e2 0;];
        Q3(1:4,1:4)=Q3_U;
        Q3(m-3:m,m-3:m)=flipud( fliplr( -Q3_U ) );
        Q3=Q3/h^2;
        d2_U=[1 -2 1;]/h^2;
        d2_l=zeros(1,m);
        d2_l(1:3)=d2_U;
        d2_r=zeros(1,m);
        d2_r(m-2:m)=fliplr(d2_U);
        D3=HI*(Q3 - e_l*d2_l + e_r*d2_r +1/2*d1_l'*d1_l -1/2*d1_r'*d1_r ) ;

        


    elseif order == 4
        H=diag(ones(m,1),0);
        H_U=[0.35809e5 / 0.100800e6 0 0 0 0 0; 0 0.13297e5 / 0.11200e5 0 0 0 0; 0 0 0.5701e4 / 0.5600e4 0 0 0; 0 0 0 0.45109e5 / 0.50400e5 0 0; 0 0 0 0 0.35191e5 / 0.33600e5 0; 0 0 0 0 0 0.33503e5 / 0.33600e5;];
        H(1:6,1:6)=H_U;
        H(m-5:m,m-5:m)=fliplr(flipud(H_U));
        H=H*h;
        HI=inv(H);
        
        % First derivative SBP operator, 1st order accurate at first 6 boundary points
        q2=-1/12;q1=8/12;
        Q=q2*(diag(ones(m-2,1),2) - diag(ones(m-2,1),-2))+q1*(diag(ones(m-1,1),1)-diag(ones(m-1,1),-1));
        Q_U = [0 0.526249e6 / 0.907200e6 -0.10819e5 / 0.777600e6 -0.50767e5 / 0.907200e6 -0.631e3 / 0.28800e5 0.91e2 / 0.7776e4; -0.526249e6 / 0.907200e6 0 0.1421209e7 / 0.2721600e7 0.16657e5 / 0.201600e6 -0.8467e4 / 0.453600e6 -0.33059e5 / 0.5443200e7; 0.10819e5 / 0.777600e6 -0.1421209e7 / 0.2721600e7 0 0.631187e6 / 0.1360800e7 0.400139e6 / 0.5443200e7 -0.8789e4 / 0.302400e6; 0.50767e5 / 0.907200e6 -0.16657e5 / 0.201600e6 -0.631187e6 / 0.1360800e7 0 0.496403e6 / 0.907200e6 -0.308533e6 / 0.5443200e7; 0.631e3 / 0.28800e5 0.8467e4 / 0.453600e6 -0.400139e6 / 0.5443200e7 -0.496403e6 / 0.907200e6 0 0.1805647e7 / 0.2721600e7; -0.91e2 / 0.7776e4 0.33059e5 / 0.5443200e7 0.8789e4 / 0.302400e6 0.308533e6 / 0.5443200e7 -0.1805647e7 / 0.2721600e7 0;];
        Q(1:6,1:6)=Q_U;
        Q(m-5:m,m-5:m)=flipud( fliplr( -Q_U ) );
        e_l=zeros(m,1);e_l(1)=1;
        e_r=zeros(m,1);e_r(m)=1;
        D1=HI*(Q-1/2*e_l*e_l'+1/2*e_r*e_r') ;

        % Second derivative, 1st order accurate at first 6 boundary points
        m2=1/12;m1=-16/12;m0=30/12;
        M=m2*(diag(ones(m-2,1),2)+diag(ones(m-2,1),-2))+m1*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1))+m0*diag(ones(m,1),0);
        M_U=[0.2386127e7 / 0.2177280e7 -0.515449e6 / 0.453600e6 -0.10781e5 / 0.777600e6 0.61567e5 / 0.1360800e7 0.6817e4 / 0.403200e6 -0.1069e4 / 0.136080e6; -0.515449e6 / 0.453600e6 0.4756039e7 / 0.2177280e7 -0.1270009e7 / 0.1360800e7 -0.3751e4 / 0.28800e5 0.3067e4 / 0.680400e6 0.119459e6 / 0.10886400e8; -0.10781e5 / 0.777600e6 -0.1270009e7 / 0.1360800e7 0.111623e6 / 0.60480e5 -0.555587e6 / 0.680400e6 -0.551339e6 / 0.5443200e7 0.8789e4 / 0.453600e6; 0.61567e5 / 0.1360800e7 -0.3751e4 / 0.28800e5 -0.555587e6 / 0.680400e6 0.1025327e7 / 0.544320e6 -0.464003e6 / 0.453600e6 0.222133e6 / 0.5443200e7; 0.6817e4 / 0.403200e6 0.3067e4 / 0.680400e6 -0.551339e6 / 0.5443200e7 -0.464003e6 / 0.453600e6 0.5074159e7 / 0.2177280e7 -0.1784047e7 / 0.1360800e7; -0.1069e4 / 0.136080e6 0.119459e6 / 0.10886400e8 0.8789e4 / 0.453600e6 0.222133e6 / 0.5443200e7 -0.1784047e7 / 0.1360800e7 0.1812749e7 / 0.725760e6;];
        M(1:6,1:6)=M_U;
        M(m-5:m,m-5:m)=flipud( fliplr( M_U ) );
        M=M/h;
        d1_U=[-0.11e2 / 0.6e1 3 -0.3e1 / 0.2e1 0.1e1 / 0.3e1;]/h;
        d1_l=zeros(1,m);
        d1_l(1:4)=d1_U;
        d1_r=zeros(1,m);
        d1_r(m-3:m)=fliplr(-d1_U);
        D2=HI*(-M-e_l*d1_l+e_r*d1_r);

        % Third derivative, 1st order accurate at first 6 boundary points
        q3=-1/8;q2=1;q1=-13/8;
        Q3=q3*(diag(ones(m-3,1),3)-diag(ones(m-3,1),-3))+q2*(diag(ones(m-2,1),2)-diag(ones(m-2,1),-2))+q1*(diag(ones(m-1,1),1)-diag(ones(m-1,1),-1));   
        Q3_U = [0 -0.88471e5 / 0.67200e5 0.58139e5 / 0.33600e5 -0.1167e4 / 0.2800e4 -0.89e2 / 0.11200e5 0.7e1 / 0.640e3; 0.88471e5 / 0.67200e5 0 -0.43723e5 / 0.16800e5 0.46783e5 / 0.33600e5 -0.191e3 / 0.3200e4 -0.1567e4 / 0.33600e5; -0.58139e5 / 0.33600e5 0.43723e5 / 0.16800e5 0 -0.4049e4 / 0.2400e4 0.29083e5 / 0.33600e5 -0.71e2 / 0.1400e4; 0.1167e4 / 0.2800e4 -0.46783e5 / 0.33600e5 0.4049e4 / 0.2400e4 0 -0.8591e4 / 0.5600e4 0.10613e5 / 0.11200e5; 0.89e2 / 0.11200e5 0.191e3 / 0.3200e4 -0.29083e5 / 0.33600e5 0.8591e4 / 0.5600e4 0 -0.108271e6 / 0.67200e5; -0.7e1 / 0.640e3 0.1567e4 / 0.33600e5 0.71e2 / 0.1400e4 -0.10613e5 / 0.11200e5 0.108271e6 / 0.67200e5 0;];
        Q3(1:6,1:6)=Q3_U;
        Q3(m-5:m,m-5:m)=flipud( fliplr( -Q3_U ) );
        Q3=Q3/h^2;
        d2_U=[2 -5 4 -1;]/h^2;
        d2_l=zeros(1,m);
        d2_l(1:4)=d2_U;
        d2_r=zeros(1,m);
        d2_r(m-3:m)=fliplr(d2_U);
        D3=HI*(Q3 - e_l*d2_l + e_r*d2_r +1/2*d1_l'*d1_l -1/2*d1_r'*d1_r ) ;

        

        
        

    elseif order == 6
        H=diag(ones(m,1),0);
        H_U=[0.318365e6 / 0.1016064e7 0 0 0 0 0 0 0; 0 0.145979e6 / 0.103680e6 0 0 0 0 0 0; 0 0 0.139177e6 / 0.241920e6 0 0 0 0 0; 0 0 0 0.964969e6 / 0.725760e6 0 0 0 0; 0 0 0 0 0.593477e6 / 0.725760e6 0 0 0; 0 0 0 0 0 0.52009e5 / 0.48384e5 0 0; 0 0 0 0 0 0 0.141893e6 / 0.145152e6 0; 0 0 0 0 0 0 0 0.1019713e7 / 0.1016064e7;];
        H(1:8,1:8)=H_U;
        H(m-7:m,m-7:m)=fliplr(flipud(H_U));
        H=H*h;
        HI=inv(H);

        % First derivative SBP operator, 1st order accurate at first 6 boundary points
        q3=1/60;q2=-3/20;q1=3/4;
        Q=q3*(diag(ones(m-3,1),3) - diag(ones(m-3,1),-3))+q2*(diag(ones(m-2,1),2) - diag(ones(m-2,1),-2))+q1*(diag(ones(m-1,1),1)-diag(ones(m-1,1),-1));
        Q_U = [0 0.1547358409e10 / 0.2421619200e10 -0.422423e6 / 0.11211200e8 -0.1002751721e10 / 0.8717829120e10 -0.15605263e8 / 0.484323840e9 0.1023865e7 / 0.24216192e8 0.291943739e9 / 0.21794572800e11 -0.24659e5 / 0.2534400e7; -0.1547358409e10 / 0.2421619200e10 0 0.23031829e8 / 0.62899200e8 0.10784027e8 / 0.34594560e8 0.2859215e7 / 0.31135104e8 -0.45982103e8 / 0.345945600e9 -0.26681e5 / 0.1182720e7 0.538846039e9 / 0.21794572800e11; 0.422423e6 / 0.11211200e8 -0.23031829e8 / 0.62899200e8 0 0.28368209e8 / 0.69189120e8 -0.9693137e7 / 0.69189120e8 0.1289363e7 / 0.17740800e8 -0.39181e5 / 0.5491200e7 -0.168647e6 / 0.24216192e8; 0.1002751721e10 / 0.8717829120e10 -0.10784027e8 / 0.34594560e8 -0.28368209e8 / 0.69189120e8 0 0.5833151e7 / 0.10644480e8 0.4353179e7 / 0.69189120e8 0.2462459e7 / 0.155675520e9 -0.215471e6 / 0.10762752e8; 0.15605263e8 / 0.484323840e9 -0.2859215e7 / 0.31135104e8 0.9693137e7 / 0.69189120e8 -0.5833151e7 / 0.10644480e8 0 0.7521509e7 / 0.13837824e8 -0.1013231e7 / 0.11531520e8 0.103152839e9 / 0.8717829120e10; -0.1023865e7 / 0.24216192e8 0.45982103e8 / 0.345945600e9 -0.1289363e7 / 0.17740800e8 -0.4353179e7 / 0.69189120e8 -0.7521509e7 / 0.13837824e8 0 0.67795697e8 / 0.98841600e8 -0.17263733e8 / 0.151351200e9; -0.291943739e9 / 0.21794572800e11 0.26681e5 / 0.1182720e7 0.39181e5 / 0.5491200e7 -0.2462459e7 / 0.155675520e9 0.1013231e7 / 0.11531520e8 -0.67795697e8 / 0.98841600e8 0 0.1769933569e10 / 0.2421619200e10; 0.24659e5 / 0.2534400e7 -0.538846039e9 / 0.21794572800e11 0.168647e6 / 0.24216192e8 0.215471e6 / 0.10762752e8 -0.103152839e9 / 0.8717829120e10 0.17263733e8 / 0.151351200e9 -0.1769933569e10 / 0.2421619200e10 0;];
        Q(1:8,1:8)=Q_U;
        Q(m-7:m,m-7:m)=flipud( fliplr( -Q_U ) );
        e_l=zeros(m,1);e_l(1)=1;
        e_r=zeros(m,1);e_r(m)=1;
        D1=HI*(Q-1/2*e_l*e_l'+1/2*e_r*e_r') ;

        % Second derivative, 1st order accurate at first 6 boundary points
        m3=-1/90;m2=3/20;m1=-3/2;m0=49/18;
        M=m3*(diag(ones(m-3,1),3)+diag(ones(m-3,1),-3))+m2*(diag(ones(m-2,1),2)+diag(ones(m-2,1),-2))+m1*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1))+m0*diag(ones(m,1),0);
        M_U=[0.4347276223e10 / 0.3736212480e10 -0.1534657609e10 / 0.1210809600e10 0.68879e5 / 0.3057600e7 0.1092927401e10 / 0.13076743680e11 0.18145423e8 / 0.968647680e9 -0.1143817e7 / 0.60540480e8 -0.355447739e9 / 0.65383718400e11 0.56081e5 / 0.16473600e8; -0.1534657609e10 / 0.1210809600e10 0.42416226217e11 / 0.18681062400e11 -0.228654119e9 / 0.345945600e9 -0.12245627e8 / 0.34594560e8 -0.2995295e7 / 0.46702656e8 0.52836503e8 / 0.691891200e9 0.119351e6 / 0.12812800e8 -0.634102039e9 / 0.65383718400e11; 0.68879e5 / 0.3057600e7 -0.228654119e9 / 0.345945600e9 0.5399287e7 / 0.4193280e7 -0.24739409e8 / 0.34594560e8 0.7878737e7 / 0.69189120e8 -0.1917829e7 / 0.31449600e8 0.39727e5 / 0.3660800e7 0.10259e5 / 0.4656960e7; 0.1092927401e10 / 0.13076743680e11 -0.12245627e8 / 0.34594560e8 -0.24739409e8 / 0.34594560e8 0.7780367599e10 / 0.3736212480e10 -0.70085363e8 / 0.69189120e8 -0.500209e6 / 0.6289920e7 -0.311543e6 / 0.17962560e8 0.278191e6 / 0.21525504e8; 0.18145423e8 / 0.968647680e9 -0.2995295e7 / 0.46702656e8 0.7878737e7 / 0.69189120e8 -0.70085363e8 / 0.69189120e8 0.7116321131e10 / 0.3736212480e10 -0.545081e6 / 0.532224e6 0.811631e6 / 0.11531520e8 -0.84101639e8 / 0.13076743680e11; -0.1143817e7 / 0.60540480e8 0.52836503e8 / 0.691891200e9 -0.1917829e7 / 0.31449600e8 -0.500209e6 / 0.6289920e7 -0.545081e6 / 0.532224e6 0.324760747e9 / 0.138378240e9 -0.65995697e8 / 0.49420800e8 0.1469203e7 / 0.13759200e8; -0.355447739e9 / 0.65383718400e11 0.119351e6 / 0.12812800e8 0.39727e5 / 0.3660800e7 -0.311543e6 / 0.17962560e8 0.811631e6 / 0.11531520e8 -0.65995697e8 / 0.49420800e8 0.48284442317e11 / 0.18681062400e11 -0.1762877569e10 / 0.1210809600e10; 0.56081e5 / 0.16473600e8 -0.634102039e9 / 0.65383718400e11 0.10259e5 / 0.4656960e7 0.278191e6 / 0.21525504e8 -0.84101639e8 / 0.13076743680e11 0.1469203e7 / 0.13759200e8 -0.1762877569e10 / 0.1210809600e10 0.10117212851e11 / 0.3736212480e10;];
        M(1:8,1:8)=M_U;
        M(m-7:m,m-7:m)=flipud( fliplr( M_U ) );
        M=M/h; 
        d1_U=[-0.25e2 / 0.12e2 4 -3 0.4e1 / 0.3e1 -0.1e1 / 0.4e1;]/h;
        d1_l=zeros(1,m);
        d1_l(1:5)=d1_U;
        d1_r=zeros(1,m);
        d1_r(m-4:m)=fliplr(-d1_U);
        D2=HI*(-M-e_l*d1_l+e_r*d1_r);

        % Third derivative, 1st order accurate at first 6 boundary points
        q4=7/240;q3=-3/10;q2=169/120;q1=-61/30;
        Q3=q4*(diag(ones(m-4,1),4)-diag(ones(m-4,1),-4))+q3*(diag(ones(m-3,1),3)-diag(ones(m-3,1),-3))+q2*(diag(ones(m-2,1),2)-diag(ones(m-2,1),-2))+q1*(diag(ones(m-1,1),1)-diag(ones(m-1,1),-1));   
        Q3_U = [0 -0.10882810591e11 / 0.5811886080e10 0.398713069e9 / 0.132088320e9 -0.1746657571e10 / 0.1162377216e10 0.56050639e8 / 0.145297152e9 -0.11473393e8 / 0.1162377216e10 -0.38062741e8 / 0.1452971520e10 0.30473e5 / 0.4392960e7; 0.10882810591e11 / 0.5811886080e10 0 -0.3720544343e10 / 0.830269440e9 0.767707019e9 / 0.207567360e9 -0.1047978301e10 / 0.830269440e9 0.1240729e7 / 0.14826240e8 0.6807397e7 / 0.55351296e8 -0.50022767e8 / 0.1452971520e10; -0.398713069e9 / 0.132088320e9 0.3720544343e10 / 0.830269440e9 0 -0.2870078009e10 / 0.830269440e9 0.74962049e8 / 0.29652480e8 -0.12944857e8 / 0.30750720e8 -0.17846623e8 / 0.103783680e9 0.68707591e8 / 0.1162377216e10; 0.1746657571e10 / 0.1162377216e10 -0.767707019e9 / 0.207567360e9 0.2870078009e10 / 0.830269440e9 0 -0.727867087e9 / 0.276756480e9 0.327603877e9 / 0.207567360e9 -0.175223717e9 / 0.830269440e9 0.1353613e7 / 0.726485760e9; -0.56050639e8 / 0.145297152e9 0.1047978301e10 / 0.830269440e9 -0.74962049e8 / 0.29652480e8 0.727867087e9 / 0.276756480e9 0 -0.1804641793e10 / 0.830269440e9 0.311038417e9 / 0.207567360e9 -0.1932566239e10 / 0.5811886080e10; 0.11473393e8 / 0.1162377216e10 -0.1240729e7 / 0.14826240e8 0.12944857e8 / 0.30750720e8 -0.327603877e9 / 0.207567360e9 0.1804641793e10 / 0.830269440e9 0 -0.1760949511e10 / 0.830269440e9 0.2105883973e10 / 0.1452971520e10; 0.38062741e8 / 0.1452971520e10 -0.6807397e7 / 0.55351296e8 0.17846623e8 / 0.103783680e9 0.175223717e9 / 0.830269440e9 -0.311038417e9 / 0.207567360e9 0.1760949511e10 / 0.830269440e9 0 -0.1081094773e10 / 0.528353280e9; -0.30473e5 / 0.4392960e7 0.50022767e8 / 0.1452971520e10 -0.68707591e8 / 0.1162377216e10 -0.1353613e7 / 0.726485760e9 0.1932566239e10 / 0.5811886080e10 -0.2105883973e10 / 0.1452971520e10 0.1081094773e10 / 0.528353280e9 0;];
        Q3(1:8,1:8)=Q3_U;
        Q3(m-7:m,m-7:m)=flipud( fliplr( -Q3_U ) );
        Q3=Q3/h^2;
        d2_U=[0.35e2 / 0.12e2 -0.26e2 / 0.3e1 0.19e2 / 0.2e1 -0.14e2 / 0.3e1 0.11e2 / 0.12e2;]/h^2;
        d2_l=zeros(1,m);
        d2_l(1:5)=d2_U;
        d2_r=zeros(1,m);
        d2_r(m-4:m)=fliplr(d2_U);
        D3=HI*(Q3 - e_l*d2_l + e_r*d2_r +1/2*d1_l'*d1_l -1/2*d1_r'*d1_r ) ;

        


        
    end



end