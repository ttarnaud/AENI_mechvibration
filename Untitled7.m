x = 0.01:0.01:10
y = 1-exp(-x)
figure
subplot(1,2,1)
for i=1:4
plot(x,y)
xlabel('lin')
subplot(1,2,2)
plot(x,y)
set(gca,'xscale','log')
xlabel('log')
mtit('1-exp(-x)')