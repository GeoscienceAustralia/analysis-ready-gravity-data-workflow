x = sort(20*rand(100,1));
v = besselj(0,x);
F = griddedInterpolant(x,v)

xq = linspace(0,20,500);
vq = F(xq);
figure(1)
plot(x,v,'ro')
hold on
plot(xq,vq,'.')
legend('Sample Points','Interpolated Values')



Fcustom = custom_griddedInterpolant({x},v)

vqcustom = Fcustom.interp(xq);

%vqcustom = Fcustom(xq);
figure(2)
plot(x,v,'ro')
hold on
plot(xq,vqcustom,'.')
legend('Sample Points','Interpolated Values')

figure(3)
plot(vqcustom-vq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x,y] = ndgrid(-5:0.8:5);
z = sin(x.^2 + y.^2) ./ (x.^2 + y.^2);
figure(4)
surf(x,y,z)

F2 = griddedInterpolant(x,y,z);

[xq,yq] = ndgrid(-5:0.1:5);
vq2 = F2(xq,yq);
figure(5)
surf(xq,yq,vq2)

Fcustom2 = custom_griddedInterpolant({x,y},z)

vqcustom2 = Fcustom2.interp(xq,yq);

figure(6)
surf(xq,yq,vqcustom2)



