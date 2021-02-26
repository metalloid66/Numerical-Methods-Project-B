%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Task 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% interval setup
interval1 = [2 6];
a1 = interval1(1);
b1 = interval1(2);
interval2 = [6 11];
a2 = interval2(1);
b2 = interval2(2);

%function set-up
f = @(x)(x*cos(x)-log(x + 2));
figure(3);
fplot(f,[2 11]);
title('f(x) = x*cos(x) – ln(x+2), [2,11]');


%                       The main program test

%start of bisection method test
[x1, itrC1] = bisection(a1,b1,f,0.0001); 
[x2, itrC2] = bisection(a2,b2,f,0.0001);
disp(['x1 = ',num2str(x1)]);
disp(['x2 = ',num2str(x2)]);
disp(['Number of iterations to get all zeros: ',num2str(itrC1+ itrC2)]);
%end of bisection method test

%start of newtons method test
[x1, itrC1] = newtons(5,f,0.0001); 
[x2, itrC2] = newtons(7,f,0.0001);
disp(['approximation point to x1 = 2']);
disp(['approximation point to x2 = 11']);
disp(['x1 = ',num2str(x1)]);
disp(['x2 = ',num2str(x2)]);
disp(['Number of iterations to get all zeros: ',num2str(itrC1+ itrC2)]);
%end of newtons method test

%                       End of the main program test



%                     Start of Newton's method Function
function [x, itr] = newtons(approx,f,acc) %takes an approximation of x and an accuracy and a function
itr = 1;
approxIni = approx;
syms x;
while 1 
    approx = approx - (f(approx)./eval(subs(diff(f(x)),approx)));
%     disp(approx);
    itr = itr+1;
    if abs(f(approx)) < acc
        x = approx;
        break
    end
end
figure(2);
fplot(f,[2 11],'r');
title('Newton''s method');
hold on 
plot(x,0,'r*'); %plotting the zero
hold on
plot(approxIni,0,'b*'); %plotting the guessed initial approximation to zero
end

%                     End of Newton's method Function

%                     Start of Bisection method Function
function [x, itr] = bisection(a,b,f,acc)
itr = 0;
fa = f(a);
while b-a > acc
    c = 0.5 * (a+b);
    fc = f(c);
    if fa * fc > 0
        a = c;
    else
        b = c;
    end
    itr = itr+1;
end
x = c;
figure(1);
fplot(f,[2 11]);
title('Bisection method');
hold on 
plot(x,0,'r*'); %pointing the zero
hold on
plot(2,0,'b*'); %pointing where the interval starts
hold on 
plot(11,0,'b*');%pointing where the interval ends
end
%                     End of Bisection method Function
