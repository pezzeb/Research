function [ out ] = PSquadmethod(fnc,A,B,modeltype)
%PSQUADMETHOD A replica of quadgk and quadl - written to understand the
%algorithm.
%   This function was developed since the functions quadgk and quadl could
%   not handle the automatic differentiation that is available to us and
%   therefore this was written. The code is inspired by the quadl and
%   quadgk implementation but not copied. Thuis function is not
%   automatic/adaptive, i.e. the intevals is fixed and no tolerance is
%   used.

%SOME LINKS:
% http://blogs.mathworks.com/cleve/2016/05/23/modernization-of-numerical-integration-from-quad-to-integral/#07a23423-d0df-4c19-a51e-3c163ad0b85f
% https://www.mathworks.com/content/dam/mathworks/mathworks-dot-com/moler/quad.pdf
% http://www.sciencedirect.com/science/article/pii/S037704270600700X

% Debugging
% global iterationVar
% iterationVar = iterationVar + 1

if(strcmpi(modeltype,'quadgk'))
    %% Gauss Knotrod motsvarande quadgk
    if(inf==B)
        I = 0;
                lu = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.90,0.9500,0.9625,0.975,0.978125,0.98125,0.984375,0.9859375,0.9875,0.98828125,0.9890625,0.98984375,0.990625,0.99140625,0.9921875,0.99296875,0.99375,1];  %lower and upper
%         lu = (0:1/650:1);
%                 lu = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
        %         lu = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.97];
        
        for i=1:length(lu)-1
            %             [It,~] = GaussKronrod(fnc,0.1*i-0.1,0.1*i,A,B);
            [It,~] = GaussKronrod(fnc,lu(i),lu(i+1),A,B);
            if(It==inf)
                I = inf;
                break;
            else
                I = I + It;
            end
        end
        %These prints are for debugging
        %         I
        %         sum(I)
        out = I;
    else
        I = 0;
        for i=1:10
            [I(i,1),Ierr(i,1)] = GaussKronrod(fnc,0.2*i-1.2,0.2*i-1.0,A,B);
        end
        %These prints are for debugging
        %         I
        %         sum(I)
        out = sum(I);
    end
elseif(strcmpi(modeltype,'quadl'))
    %% Gauss Kronrod motsvarande quadl - MEN TOLERANSEN �R INTE DEN SAMMA...DEN �R H�RDKODAD H�R
    I1 = GaussLobbato(fnc,A,B,4);
    I2 = LobbatoKronrod(fnc,A,B);
    
    if(abs(I1-I2)>=0.007705959723185)
        I = 0;
        
        [I1,I2] = LobbatoKronrodhelp(fnc,A,B,0.007705959723185);
    else
        %nothing
    end
    out = I2;
end
end

function [I1,I2] = LobbatoKronrodhelp(fnc,a,b,tol)
I1 = GaussLobbato(fnc,a,b,4);
I2 = LobbatoKronrod(fnc,a,b);
if(abs(I1-I2)<tol)
    return;
else
    c = (a + b)/2;
    h = (b - a)/2;
    alpha = sqrt(2/3);
    beta = 1/sqrt(5);
    x = [c-alpha*h c-beta*h c c+beta*h c+alpha*h];
    x = [a x b];
    I=0;
    for k=1:6
        [~,I2] = LobbatoKronrodhelp(fnc,x(k),x(k+1),tol);
        I = I + I2;
    end
    I2 = I;
end
end
function out = midpointMethod(fnc,a,b)
out = (b - a)*fnc((a + b)/2);
end
function out = trapezoidMethod(fnc,a,b)
out = (b - a)/2*(fnc(a)+fnc(b));
end
function out = SimpsonMethod(fnc,a,b)
out = (b - a)/6*(fnc(a) + 4*fnc((a+b)/2) + fnc(b));
end
function out = GaussLegend(fnc,a,b)

fnc_11 = @(t) fnc( (b-a)/2*t + (a + b)/2 );
%three points
out = 8/9*fnc_11(0) + 5/9*( fnc_11(sqrt(3/5)) + fnc_11(-sqrt(3/5) ) );

end
function out = GaussLobbato(fnc,a,b,n)

fnc_11 = @(t) fnc( (b-a)/2*t + (a + b)/2 );
%three points

if(3==n)
    x = [0;-1;1];
    w = [4/3;1/3;1/3];
elseif(4==n)
    x = [sqrt(1/5);-sqrt(1/5);-1;1];
    w = [5/6;5/6;1/6;1/6];
end
out = (b-a)/2*w'*fnc_11(x);

end
function out = LobbatoKronrod(fnc,a,b)

c = (a + b)/2;
h = (b - a)/2;
alpha = sqrt(2/3);
beta = 1/sqrt(5);
x = [c-alpha*h c-beta*h c c+beta*h c+alpha*h];
x = [a x b];
y = feval(fnc,x);
out = (h/1470)*[77 432 625 672 625 432 77]*y.';

end
function [Iknot,Ierr] = GaussKronrod(fnc,a,b,A,B)
%% REFERENCES AND EXPLANATION OF THIS FUNCTION
% This function is based on the quadgk implementation but this version is
% static and is thus not an adaptive/automatic implementation. Thus, a good
% method of understanding this method is by studying the quadgk
% implementation. The quadgk implementation can also handle waypoints which
% this algorithm does not handle. Another good (first) reference is the
% book: Scientific Computing an Introductory Survey 2nd ed.

% The quadgk and the function integral have the same implementation:
% http://blogs.mathworks.com/cleve/2016/05/23/modernization-of-numerical-integration-from-quad-to-integral/#07a23423-d0df-4c19-a51e-3c163ad0b85f

% There are four cases, i.e. for different integrationsintervals. [a,b],
% [a,inf], [-inf,b] and [-inf,inf]. But here only two of the are
% implemented. but see quadgk for the rest.

% An important part of the algorithm is the weights for Gauss and Knotrod.
% These can be found in a number of reference. One is: https://www.advanpix.com/2011/11/07/gauss-kronrod-quadrature-nodes-weights/
% an other option is to instead look into quadgk but there is a bit harder
% to find but possible.
%
% Regarding the weights. The weights for Kronrod can be reused for the
% gaussian or perhaps more correct. The Kronrod weights is the gauss
% weights plus some more. And the same holds for the points. So instead of
% calculating the gauss weights we can calcualte the difference between
% kronrod and gauss directly.

% A transformation is used for weaken problems in the end points and also
% to remove some other problems. Read a bit about it in the book: Methods
% of Numerical Integration page 441. Another good article, Ref: L.F.
% Shampine, "Vectorized Adaptive Quadrature in Matlab", Journal of
% Computational and Applied Mathematics 211, 2008, pp.131-140. This article
% is also referenced in the quadgk implementation.


%% Code start

%Faktorn: (b-a)/2 �r det som i quadgk rad 257 �r halfh.
if(B~=inf)
    % fnc_11 = @(t) (b-a)/2*fnc( (b-a)/2*t + (a + b)/2 ); %Ps version - Denna version anv�nder inte en transform
    
    fnc_11 = @(y) 3/4.*(B-A).*(1 - y.^2).* fnc( 1/4.*(B-A).*y.*(3-y.^2) + 1/2.*(B+A) ).*(b-a)/2;
else
    fnc_11 = @(y) 2.*y./(1 - y)./((1 - y).^2).*fnc((y./(1 - y)).^2 + A).*(b-a)/2;
end

%Points (or more correct: Abscissas)
xk_nod = [ -9.9145537112081263920685469752632852e-01;...
    -9.4910791234275852452618968404785126e-01;...
    -8.6486442335976907278971278864092620e-01;...
    -7.4153118559939443986386477328078841e-01;...
    -5.8608723546769113029414483825872960e-01;...
    -4.0584515137739716690660641207696146e-01;...
    -2.0778495500789846760068940377324491e-01;...
    1.8367099231598242312011508394097589e-40;...
    2.0778495500789846760068940377324491e-01;...
    4.0584515137739716690660641207696146e-01;...
    5.8608723546769113029414483825872960e-01;...
    7.4153118559939443986386477328078841e-01;...
    8.6486442335976907278971278864092620e-01;...
    9.4910791234275852452618968404785126e-01;...
    9.9145537112081263920685469752632852e-01];

xk = (b-a)/2*xk_nod + (a+b)/2;

%Weights
wk = [2.2935322010529224963732008058969592e-02;...
    6.3092092629978553290700663189204287e-02;...
    1.0479001032225018383987632254151802e-01;...
    1.4065325971552591874518959051023792e-01;...
    1.6900472663926790282658342659855028e-01;...
    1.9035057806478540991325640242101368e-01;...
    2.0443294007529889241416199923464908e-01;...
    2.0948214108472782801299917489171426e-01;...
    2.0443294007529889241416199923464908e-01;...
    1.9035057806478540991325640242101368e-01;...
    1.6900472663926790282658342659855028e-01;...
    1.4065325971552591874518959051023792e-01;...
    1.0479001032225018383987632254151802e-01;...
    6.3092092629978553290700663189204287e-02;...
    2.2935322010529224963732008058969592e-02];

wg = [ 1.2948496616886969327061143267908202e-01;...
    2.7970539148927666790146777142377958e-01;...
    3.8183005050511894495036977548897513e-01;...
    4.1795918367346938775510204081632653e-01;...
    3.8183005050511894495036977548897513e-01;...
    2.7970539148927666790146777142377958e-01;...
    1.2948496616886969327061143267908202e-01];

%Calculate the integral and the integral error
ftmp  = fnc_11(xk);
if(not(isa(ftmp,'forADm')) && not(isa(ftmp,'revADm')))
    ftmp(isnan(ftmp))=0;
end
if(any(ftmp==inf)||any(ftmp==-inf))
    Iknot = inf;
    Ierr  = inf;
else
    Iknot = wk'*ftmp;
    
    ewk = wk;
    ewk(2:2:14) = ewk(2:2:14)-wg;
    
    Ierr = ewk'*fnc_11(xk);
end
end
%% lglnodes - DENNA FUNKTION �R TAGEN FR�N INTERNET: https://se.mathworks.com/matlabcentral/fileexchange/4775-legende-gauss-lobatto-nodes-and-weights?focused=5055264&tab=function

%OBSERVERA att jag har gjort en f�r�ndringen som f�reslogs av mathworks
%kommentarsf�lt. Dvs x=-cos(pi*(0:N)/N)'; har ett minustecken, vilket inte
%finns i orginal implementeringen.

function [x,w,P]=lglnodes(N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% lglnodes.m
%
% Computes the Legendre-Gauss-Lobatto nodes, weights and the LGL Vandermonde
% matrix. The LGL nodes are the zeros of (1-x^2)*P'_N(x). Useful for numerical
% integration and spectral methods.
%
% Reference on LGL nodes and weights:
%   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
%   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
%
% Written by Greg von Winckel - 04/17/2004
% Contact: gregvw@chtm.unm.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Truncation + 1
N1=N+1;

% Use the Chebyshev-Gauss-Lobatto nodes as the first guess
x=-cos(pi*(0:N)/N)';

% The Legendre Vandermonde Matrix
P=zeros(N1,N1);

% Compute P_(N) using the recursion relation
% Compute its first and second derivatives and
% update x using the Newton-Raphson method.

xold=2;

while max(abs(x-xold))>eps
    
    xold=x;
    
    P(:,1)=1;    P(:,2)=x;
    
    for k=2:N
        P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
    end
    
    x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
    
end

w=2./(N*N1*P(:,N1).^2);
end