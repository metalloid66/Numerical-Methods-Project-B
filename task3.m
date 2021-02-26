%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Task 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function set-up
syms x;
ff = @(x)(2*x.^4+3*x.^3-5*x.^2+2*x+1);
firstDiff = @(x)(8*x.^3 + 9*x.^2 - 10*x + 2);
secondDiff = @(x)(24*x.^2 + 18*x - 10);
n = polynomialDegree(ff(x));

%                       The main program test

[x1Lag,itrC1] = laguerre(-1,ff,n,firstDiff,secondDiff,0.000001);
[x2Lag,itrC2] = laguerre(-3,ff,n,firstDiff,secondDiff,0.000001);
[x3Lag,itrC3] = laguerre(1,ff,n,firstDiff,secondDiff,0.000001);
 x4Lag = conj(x3Lag);
disp(['Laguerre method']);
disp(['approximation point to x1 = -1']);
disp(['approximation point to x2 = -3']);
disp(['approximation point to complex x3 and x4 = 1']);
disp(['x1 = ',num2str(x1Lag)]);
disp(['x2 = ',num2str(x2Lag)]);
disp(['x3 = ',num2str(x3Lag)]);
disp(['x4 = ',num2str(x4Lag)]);
disp(['Number of iterations to get only real roots: ',num2str(itrC1+ itrC2)]);
disp(['Number of iterations to get all roots: ',num2str(itrC1+ itrC2 +itrC3)]);

%                       The end program test


%                   Start of laguerre method function

function [x,itr] = laguerre(xk,ff,n,firstDiff,secondDiff,prec)
zmin = 0;
itr = 0;
while 1 
    xk=xk+zmin;
    a =secondDiff(xk)/2;
    b =firstDiff(xk);
    c =ff(xk);
    if abs(b+sqrt(b^2-4*a*c))>=abs(b-sqrt(b^2-4*a*c))
        zmin = -n*c/(b+sqrt((n-1)*((n-1)*(b^2-4*a*c))));
    else
        zmin = -n*c/(b-sqrt((n-1)*((n-1)*(b^2-4*a*c))));
    end   
    x = xk;
    itr = itr+1;
    if abs(ff(xk))<prec
        break;
    end
end
figure(4);
fplot(ff,[-4 3],'r');
title('Laguerre method');
hold on
plot(x,0,'r*');
hold on
plot(-1,0,'b*');
hold on 
plot(-3,0,'g*');
hold on 
plot(1,0,'k*');
end
%                   End of laguerre method function