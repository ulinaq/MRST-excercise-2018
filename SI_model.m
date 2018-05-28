clear all;

i0 = 0.1;
s0 = 1 - i0;
Beta = 0.5;

% ODE
ODEfun = @(x) Beta*x*(1-x);
ODE2 = @(x) -Beta*x*(1-x);
% Solution
Sol = @(x) 1/(1+(1/i0-1)*exp(-1*Beta*x));

% Euler explicit

h = [5 1 0.2 0.04 0.008 0.0016];
tp = 10;
for i=1:length(h)
    in = i0;
    it = in;
    sn = s0;
    t = 0;
    j = 1;
    er = [];
while it<1-1e-4 & it(end)>0
in(end+1) = in(end) + Beta*in(end)*(1-in(end))*h(i);
sn(end+1) = sn(end) - Beta*sn(end)*(1-sn(end))*h(i) ;
t = t+h(i);
er(j) = abs(Sol(t)-in(end));
j = j+1;
it = in(end);
end
figure()
plot(linspace(0,t,length(sn)),sn)
figure()
plot(er)
errEx{i} = er;
end
% Error
for i=1:length(h)
    a(i)=errEx{i}(5/h(i));
end
for i=1:length(h)-1
Conv(i)=log(errEx{i}(5/h(i))/errEx{i+1}(5/h(i+1)))./log(h(i)/h(i+1));
end
% Euler implicit

for i=1:length(h)
    in = i0;
    it = in;
    sn = s0;
    t = 0;
    j = 1;
    er = [];
while it<1-1e-4 & it(end)>0
    t = t+h(i);
    in(end+1) = fzero(@(y) y - h(i)*ODEfun(y) - in(end), in(end));
%     sn(end+1) = fzero(@(y) y - h(i)*ODEfun(y) - sn(end), sn(end));
    er(j) = abs(Sol(t)-in(end));
    j = j+1;
    it = in(end);
end
figure()
plot(linspace(0,t,length(in)),in)
figure()
plot(er)

errim{i} = er;
end
% Error
for i=1:length(h)-1
    ConvIm(i)=log(errim{i}(5/h(i))/errim{i+1}(5/h(i+1)))./log(h(i)/h(i+1));
end
for i=1:length(h)
    b(i)=errim{i}(5/h(i));
end
% Error
