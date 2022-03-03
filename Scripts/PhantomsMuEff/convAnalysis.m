res = ao.result;

figure();
ax(1) = subplot(1,2,1)
plot(res.export.signal(1:4016));
title('Signal')
ax(2) = subplot(1,2,2)
plot(res.export.deMul(1:8032));
title('DeMul')

figure();
plot(res.phi);

phiConvKernel = res.phi;


% resNoHad = ao.result;
