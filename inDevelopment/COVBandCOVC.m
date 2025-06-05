% this is how to test COVA

haversineDistance=(0:pi/180/2: pi);

for i=1:1
cov0 =   COVA(300, num2str(i), haversineDistance,0, 0);
cov500 = COVA(300, num2str(i), haversineDistance,500*10^3, 500*10^3);

figure
%plot(rad2deg(haversineDistance),cov0,'r-')
hold on 
plot(rad2deg(haversineDistance),cov500,'r*')
hold on 
plot(rad2deg(haversineDistance),0*haversineDistance,'b.')
title(['covariance model ',num2str(i)])
end































