function TestNumInt

Trap(1) = 0;
Simp(1) = 0;
GL2(1) = 0;
GL3(1) = 0;

for n = 6 : 6 : 1000

    X = 0:(pi/n):pi;
    Y = exp(X).*(cos(X));
    h = pi/n;

    n = size(X,2);

    %% Trapezoidal Rule

    T = 0;
    
    for i = 1 : n
        if (i == 1 || i == n)
            T = T + Y(i)/2;
        else
            T = T + Y(i);
        end
    end
    
    T = T*h;

    Trap(round(n/6)) = T;

    %% Simpson's Rule

    S = 0;

    for i = 2:n-1
        if (mod(i,2) == 0)
            S = S + 4*Y(i);
        else
            S = S + 2*Y(i);
        end
    end

    S = S + Y(1) + Y(n);
    S = S*h/3;
    
    Simp(round(n/6)) = S;

    %% Gauss-Legendre n = 2

    n = n - 1;
    G2 = 0;
     
    l = n/2;
    hp = pi/l;
    h = hp/2;
    
    for i = 0 : l - 1 
        xi = i * hp;
        xi1 = (i+1) * hp;
        left = (xi + xi1 - (2*h)/sqrt(3))/2;
        right = (xi + xi1 + (2*h)/sqrt(3))/2;
        result = exp(left)*cos(left) + exp(right)*cos(right);
        G2 = G2 + result;
    end
    
    G2 = h*G2;
   
    GL2(round((n+1)/6)) = G2;

    %% Gauss-Legendre n = 3

    G3 = 0;

    l = n/3;
    hp = pi/l;
    h = hp/3;
    
    for i = 0 : l - 1   
        xi = i * hp;
        xi1 = (i+1) * hp;
        left = (xi + xi1 - (3*h)*sqrt(3/5))/2;
        middle = (xi + xi1)/2;
        right = (xi + xi1 + (3*h)*sqrt(3/5))/2;
        result = (5/9)*exp(left)*cos(left) + (8/9)*exp(middle)*cos(middle) + (5/9)*exp(right)*cos(right);
        G3 = G3 + result;
    end
    
    G3 = 3*h*G3/2;
    
    GL3(round((n+1)/6)) = G3;

end

TV = (1/2)*(-1-exp(pi));

Trap = abs(Trap - TV);
Simp = abs(TV - Simp);
GL2 = abs(TV - GL2);
GL3 = abs(TV - GL3);

R = 6:6:1000;

loglog(R,Trap,  'r', 'linewidth', 2); grid on; hold on;
loglog(R,Simp, 'g', 'linewidth', 2); hold on;
loglog(R,GL2, 'b', 'linewidth', 2); hold on;
loglog(R,GL3, 'k', 'linewidth', 2);

Q = GL3(GL3 > 10e-14);

[a0 a1] = LeastSquare(log(R), log(Trap));
[b0 b1] = LeastSquare(log(R), log(Simp));
[c0 c1] = LeastSquare(log(R), log(GL2));
[d0 d1] = LeastSquare(log(R(1:size(Q,2))), log(Q));


legend('show')
legend('Tn', 'Sn', 'GL2n', 'GL3n')
xlim([6 1000])
set(gca, 'XTick', [6 10 100 1000])
set(gca, 'YTick', [10^-15 10^-12 10^-8 10^-4 1])
title('f(x) = e^xcos(x)', 'FontSize', 20);
xlabel('n', 'FontSize', 14);
ylabel('Error', 'FontSize', 14);
x1 = 1000;
y1 = Trap(166); y2 = Simp(166); y3 = GL2(166); y4 = GL3(166);
t1 = num2str(a1); t2 = num2str(b1); t3 = num2str(c1); t4 = num2str(d1);
text(x1,y1,t1, 'Color', 'r'); text(x1,y2,t2, 'Color', 'g'); 
text(x1,y3,t3, 'Color', 'b'); text(x1,y4,t4, 'Color', 'k');

end

