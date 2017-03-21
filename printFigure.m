
n = 10;
k = 40;
x = 1:n;
yf = zeros(1,n);
yb = zeros(1,n);
noisestep = 0.01;
for i = 1:n;
    noise = i*noisestep;
    for j = 1:k;
        [score, ~] = stabilityMeasure(noise, 'forward');
        yf(i) = yf(i) + score;
        [score, ~] = stabilityMeasure(noise, 'both');
        yb(i) = yb(i) + score;
    end
    yf(i) = yf(i)/k;
    yb(i) = yb(i)/k;
end
figure('Visible','on');
plot( x, log(yf), x, log(yb));
legend('Fitting Eigenvector','Avoiding Fit of Eigenvector');
xlabel('Noise Level');
ylabel('Score (Log Scale');
