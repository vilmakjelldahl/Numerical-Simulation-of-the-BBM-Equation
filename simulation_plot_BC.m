 % plot
% uu = U(1:mx,1);
% 
% hold on
% ylim([-20 37])
% title('')
% xlabel('')
% ylabel('')
% p = plot(xvec, uu, 'b');
% p.XDataSource = 'xvec';
% p.YDataSource = 'uu';
% 
% for i=1:mt*n_blocks-1
%     if mod(i,mt*50) == 0
%         uu = U(i*mx+1:(i+1)*mx,1);
% 
%         refreshdata
%         drawnow
%         % pause(0.5)
%     end
% end





%PROJEKTION
hold on
subplot(2,2,1)
plot(xvec, U(:,46))
ylim([-8 10])
title(['Numerical solution at t = ' num2str(0.01*45)])

subplot(2,2,2)
plot(xvec, U(:,51))
ylim([-8 10])
title(['Numerical solution at t = ' num2str(0.01*50)])

subplot(2,2,3)
plot(xvec, U(:,56))
ylim([-8 10])
title(['Numerical solution at t = ' num2str(0.01*55)])

subplot(2,2,4)
plot(xvec, U(:,61))
ylim([-8 10])
title(['Numerical solution at t = ' num2str(0.01*60)])







%SAT
% hold on
% subplot(2,2,1)
% plot(xvec, U(1:mx, 2000))
% ylim([-6 11])
% title(['Numerical solution at t = ' num2str(ht*2000)])
% 
% subplot(2,2,2)
% plot(xvec, U(1:mx, 2250))
% ylim([-6 11])
% title(['Numerical solution at t = ' num2str(ht*2250)])
% 
% subplot(2,2,3)
% plot(xvec, U(1:mx, 2500))
% ylim([-6 11])
% title(['Numerical solution at t = ' num2str(ht*2500)])
% 
% subplot(2,2,4)
% plot(xvec, U(1:mx, 2750))
% ylim([-6 11])
% title(['Numerical solution at t = ' num2str(ht*2750)])

