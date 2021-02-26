%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Task 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function set-up
ff = @(x)(2*x.^4+3*x.^3-5*x.^2+2*x+1);
firstDiff = @(x)(8*x.^3 + 9*x.^2 - 10*x + 2); % Derivatives to be used in newton's method
secondDiff = @(x)(24*x.^2 + 18*x - 10);


%                       The main program test

% x =-3:0.01:3;     % setting up the plot of the polynomial
% f =ff(x);      
% figure(1);
% plot(x,f);grid on;  
% axis([-5 5 -5 5])
% title('f(x) = 2x^4 + 3x^3 - 5x^2 + 2x + 1');


%start of the muller MM1 method test
[x1MM1,itrC1] = MM1(-1,-2,-3,ff,0.000001);
[x2MM1,itrC2] = MM1(1.1,-0.5,0,ff,0.000001); 
[x3MM1,itrC3] = MM1(0.1,1,2,ff,0.000001);
x4MM1 = conj(x3MM1); 

disp(['MM1 method']);
disp(['approximation points to x1 = -1,-2,-3']);
disp(['approximation points to x2 = -1,-0.5,0']);
disp(['approximation points to complex x3 and x4 = 0.1,1,2']);
disp(['x1 = ',num2str(x1MM1)]);
disp(['x2 = ',num2str(x2MM1)]);
disp(['x3 = ',num2str(x3MM1)]);
disp(['x4 = ',num2str(x4MM1)]);
disp(['Number of iterations to get all zeros: ',num2str(itrC1+ itrC2 +itrC3)]);

%start of the muller MM2 method test
[x1MM2,itrC1] = MM2(-1,ff,firstDiff,secondDiff,0.000001);
[x2MM2,itrC2] = MM2(-3,ff,firstDiff,secondDiff,0.000001);
[x3MM2,itrC3] = MM2(1,ff,firstDiff,secondDiff,0.000001);
 x4MM2 = conj(x3MM2);
disp(['MM2 method']);
disp(['approximation point to x1 = -1']);
disp(['approximation point to x2 = -3']);
disp(['approximation point to complex x3 and x4 = 1']);
disp(['x1 = ',num2str(x1MM2)]);
disp(['x2 = ',num2str(x2MM2)]);
disp(['x3 = ',num2str(x3MM2)]);
disp(['x4 = ',num2str(x4MM2)]);
 disp(['Number of iterations to get only real roots: ',num2str(itrC1+ itrC2)]);
 disp(['Number of iterations to get all roots: ',num2str(itrC1+ itrC2 +itrC3)]);

%start of the muller newtons method test
[x1newt, itrC1] = newtons(-1,ff,0.000001); 
[x2newt, itrC2] = newtons(-3,ff,0.000001);
disp(['Newton''s method']);
disp(['approximation point to x1 = -1']);
disp(['approximation point to x2 = -3']);
disp(['x1 = ',num2str(x1newt)]);
disp(['x2 = ',num2str(x2newt)]);
disp(['Number of iterations to get all real roots: ',num2str(itrC1+ itrC2)]);

%                       End of the main program test


%                     Start of Muller MM1 method Function

function [x, itr] = MM1(farPoint,lastAprox,curAprox,ff,prec)
syms X Y;
itr = 0;
while 1
    z0 = farPoint - curAprox;
    z1 = lastAprox - curAprox;
    equ1 = z0.^2*X + z0*Y == ff(farPoint) - ff(curAprox);
    equ2 = z1.^2*X + z1*Y == ff(lastAprox) - ff(curAprox);
    [A,B] = equationsToMatrix([equ1,equ2],[X,Y]);
    resEqus = linsolve(A,B);
    a = double(resEqus(1,1));
    b = double(resEqus(2,1));
    c = ff(curAprox);
    if abs(b+sqrt(b^2-4*a*c))>=abs(b-sqrt(b^2-4*a*c))
        zmin = -2*c/(b+sqrt(b^2-4*a*c));
    else
        zmin = -2*c/(b-sqrt(b^2-4*a*c));
    end
    farPoint = lastAprox;%x1
    lastAprox = curAprox;%x2
    curAprox = curAprox + zmin; %x3
    x = curAprox;
    itr = itr+1;
    if abs(ff(curAprox))<prec
        break;
    end
end
figure(1);
fplot(ff,[-4 3],'r');
title('MM1 method');
hold on 
plot(x,0,'r*');
hold on
plot(-1,0,'b*');
hold on 
plot(-2,0,'b*');    % to plot x1
hold on 
plot(-3,0,'b*');

hold on
plot(-0.5,0,'g*');
hold on 
plot(0,0,'g*');         % to plot x2
hold on 
plot(1.1,0,'g*');

hold on
plot(0.1,0,'k*');
hold on 
plot(1,0,'k*');         % to plot x3
hold on 
plot(2,0,'k*');

end
%                     End of Muller MM1 method Function


%                     Start of Muller MM2 method Function
function [x,itr] = MM2(x0,ff,firstDiff,secondDiff,prec)
zmin =0;
xk = x0;
itr = 0;
while 1 
    xk=xk+zmin;
    a =secondDiff(xk)/2;
    b =firstDiff(xk);
    c =ff(xk);
    
    if abs(b+sqrt(b^2-4*a*c))>=abs(b-sqrt(b^2-4*a*c))
        zmin = -2*c/(b+sqrt(b^2-4*a*c));
    else
         zmin = -2*c/(b-sqrt(b^2-4*a*c));
    end
        x =xk;
        itr = itr + 1;
    if abs(ff(xk))<prec
        break;
    end
end

figure(2);
fplot(ff,[-4 3],'r');
title('MM2 method');
hold on
plot(x,0,'r*');
hold on
plot(-1,0,'b*');
hold on 
plot(-3,0,'g*'); 
hold on 
plot(1,0,'k*');
end

%                     End of Muller MM2 method Function


%                     Start of Newtons method Function

function [x, itr] = newtons(approx,f,acc) %takes an approximation of x and an accuracy
itr = 0;
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
figure(3);
fplot(f,[-4 1],'r');
title('Newton''s method');
hold on 
plot(x,0,'r*');
hold on
plot(-1,0,'b*');
hold on 
plot(-3,0,'b*');
end

%                     End of Newtons method Function