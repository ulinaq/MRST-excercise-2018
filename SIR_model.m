clear all;

i0 = 0.1;
s0 = 1 - i0;
Beta = 1;
sigma = 2;

% ODE
ODEfun = @(i,s) -1*Beta*i*s;
ODE2 = @(i,s)  Beta*i*s-Beta/sigma*i;
ODE3 = @(i) Beta/sigma*i;

% Euler explicit

h = 0.05;
tp = 10;

    in(1) = i0;
    sn(1) = s0;
    t = 0;
for j=1:700
in(j+1) = in(j) + (Beta*in(j)*sn(j)-Beta/sigma*in(j))*h;
sn(j+1) = sn(j) - h*Beta*in(j)*sn(j);
end


% Euler implicit

    ini(1) = i0;
    sni(1) = s0;
    t = 0;
for j=1:700
ini(j+1) = fzero(@(y) y - h*ODE2(y,sni(j)) - ini(j), ini(j));
sni(j+1) = fzero(@(y) y - h*ODEfun(ini(j),y) - sni(j), sni(j));
end

x=0:h:h*(length(sn)-1);
xi=0:h:h*(length(sni)-1);

figure()
plot(x,sn);
figure()
subplot(1,2,1)
plot(xi,ini,xi,sni);
subplot(1,2,2)
plot(xi,ini2,xi,sni2);
ini2=ini;sni2=sni;
in+sn-1/sigma*log(sn);
imax = i0+s0- 1/sigma - 1/sigma*(log(sigma*s0))
smax=@(x) i0 +s0 -x +1/sigma*log(x/s0);
