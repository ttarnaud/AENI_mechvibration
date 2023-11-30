% with this script we test the creation of a artificial PSD curve as input
% as varying dipole moments
x = logspace(-3,10,5000);
wn = 10;
theta = 0.2
w2 = 10^3;
alpha1 = rand()+2
alpha2 = 2*rand()+2
y = @(x) (1+x/(theta*wn))./((1+(x/wn).^alpha1).*(1+(x/w2).^alpha2))
figure(100)
subplot(2,1,1)
plot(x,10*log10(y(x)))
hold on
scatter([theta*wn,wn,w2],10*log10(y([theta*wn,wn,w2])))
[maxy,idxy] = max(10*log10(y(x)));
scatter(x(idxy),maxy)
hold off
set(gca,'xscale','log')
ylim([-200,10])
xlim([1,10^10])

x = logspace(-3,10,5000);
%x=abs_freqs;
a=2
b = 0.5;
noise1 = 10.^(a*(rand(1,length(x))-0.5))
noise2 = @(x) 10.^(b*(randn(1,length(x))))
noise3 = @(x) lognrnd(0,0.5,1,length(x))
y = @(x) (1+x/(theta*wn))./((1+(x/wn).^alpha1).*(1+(x/w2).^alpha2))
y2 = @(x) y(x).*noise2(x)
figure(100)
subplot(2,1,2)
plot(x,10*log10(y(x).*noise1))
hold on
plot(x,10*log10(y(x).*noise2(x)))
plot(x,10*log10(y(x).*noise3(x)));
%plot(x,10*log10(y(x)).*(1+(2*rand(1,length(x)))-1))
plot(x,10*log10(y(x)))
hold off

%hold off
set(gca,'xscale','log')
ylim([-200,10])
xlim([1,10^10])

%% find max
thetas = 0.3:0.1:1;
alphas = 2:0.1:3;
x = logspace(-3,10,5000);
wn = 14;

w2 = 10^3;

alpha2 = 2*rand()+2
for ith = 1:length(thetas)
    for ia = 1:length(alphas)
theta = thetas(ith);
alpha1 = alphas(ia);
y = @(x) (1+x/(theta*wn))./((1+(x/wn).^alpha1).*(1+(x/w2).^alpha2));
yval = y(x);
[maxy,idxmaxy] = max(yval);
maxx = x(idxmaxy);



THETA(ith,ia) = theta;
ALPHA(ith,ia) = alpha1;
MAXY(ith,ia) = maxy;
MAXX(ith,ia) = maxx;
    end
end

figure(200)
subplot(1,2,1)
surf(THETA,ALPHA,MAXY)
xlabel('theta')
ylabel('alpha')
subplot(1,2,2)
surf(THETA,ALPHA,MAXX)
xlabel('theta')
ylabel('alpha')