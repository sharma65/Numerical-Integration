function [ b0, b1 ] = LeastSquare(X, Y)

n = size(X,2)

xi = X;
yi = Y;

i = 1;

b0 = 0;
b1 = 0;
sumxi = sum(xi);
sumyi = sum(yi);
sumxiyi = sum(xi.*yi);
sumx2i = sum(xi.*xi);

b1 = (n*sumxiyi - sumxi*sumyi)/(n*sumx2i - sumxi*sumxi);
b0 = (sumyi/n) - b1*(sumxi/n);
