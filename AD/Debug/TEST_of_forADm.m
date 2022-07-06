%% INIT AND OTHER==========================================================
%% INIT-LOADING
clear

addpath('.\class');
% addpath('..\..\function');
addpath('.\test_function');

%A USER INTERFACE WHERE A {} IS CREATE THAT STORES ALL THE FUNCTION TO BE
%TESTED - EACH ELEMENT IN {} IS IN ITSELF A {} THAT INCLUDES A FUNCTION
%VALUE AND A VALUE AND ALSO A VECTOR OF DERIVATIVES

%OUTPUT FILE **************************************************************
% fid1 = fopen('C:\ptemp\DropboxPortableAHK\Dropbox\UNI\Exjobb\Work\implementation\AD\Debug\test_function\output\output_from_TEST_of_revADm.txt','w');
% fid2 = fopen('C:\ptemp\DropboxPortableAHK\Dropbox\UNI\Exjobb\Work\implementation\AD\Debug\test_function\output\output_from_TEST_of_revADm_matrix_operations.txt','w');
%fid1 = fopen('C:\Users\Pontus\Dropbox\UNI\Exjobb\Work\implementation\AD\Debug\test_function\output\output_from_TEST_of_revADm.txt','w');
%fid2 = fopen('C:\Users\Pontus\Dropbox\UNI\Exjobb\Work\implementation\AD\Debug\test_function\output\output_from_TEST_of_revADm_matrix_operations.txt','w');
% fid1 = fopen('C:\Users\Söderbäck\SkyDrive\Dropbox\UNI\Exjobb\Work\implementation\AD\Debug\test_function\output\output_from_TEST_of_revADm.txt','w');
% fid2 = fopen('C:\Users\Söderbäck\SkyDrive\Dropbox\UNI\Exjobb\Work\implementation\AD\Debug\test_function\output\output_from_TEST_of_revADm_matrix_operations.txt','w');
fid1 = fopen('.\outputFor\output_from_TEST_of_forADm.txt','w');
fid2 = fopen('.\outputFor\output_from_TEST_of_forADm_matrix_operations.txt','w');
%**************************************************************************
%% INPUT VARIABLES - OPERATIONS

xs_val = 3;
ys_val = 2;
zs_val = 3;
ws_val = 6;
qs_val = 12;

Bc5_val =  [1;2;3;4;5];
Cc3_val = [1;2;3];
Cr3_val = [1 2 3];
Dc3_val = [6;3;1];
Er4_val = [1 2 3 4];
Fc5_val = [3;5;6;12;2];
Gr4_val = [9 2 1 3];

xs = forADm(xs_val,[1 0 0 0 0]);
ys = forADm(ys_val,[0 1 0 0 0]);
zs = forADm(zs_val,[0 0 1 0 0]);
ws = forADm(ws_val,[0 0 0 1 0]);
qs = forADm(qs_val,[0 0 0 0 1]);

Bc5 = forADm(Bc5_val ,eye(5));
Cc3 = forADm(Cc3_val ,eye(3));
Cr3 = forADm(Cr3_val,eye(3));
Dc3  = forADm(Dc3_val ,eye(3));
Er4  = forADm(Er4_val ,[eye(4);zeros(1,4)]);
Fc5  = forADm(Fc5_val ,eye(5));
Gr4  = forADm(Gr4_val ,[eye(4);zeros(1,4)]);
%% MATRIX ARITHMETIC OPERATION - FIRST TEST SERIES=========================
%% SCALAR TEST
%function to test                   Save in as one element                 func     in  corr val          corr derivative
test_fun = @(x) x + 2;                  element1{1}                        = {test_fun,{xs},test_fun(xs_val)          ,{[1 0 0 0 0]}};
test_fun = @(x) x + 6.42;               element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[1 0 0 0 0]}};
test_fun = @(x) x - 6.42;               element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[1 0 0 0 0]}};
test_fun = @(x) 6.9 - x;                element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[-1 0 0 0 0]}};
test_fun = @(x) 6.9 + x;                element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[1 0 0 0 0]}};

test_fun = @(x) -x;                     element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[-1 0 0 0 0]}};
test_fun = @(x) 5.5*x;                  element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[5.5 0 0 0 0]}};
test_fun = @(x) x*3.4;                  element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[3.4 0 0 0 0]}};
test_fun = @(x) 7.2.*x;                 element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[7.2 0 0 0 0]}};
test_fun = @(x) x.*1.5;                 element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[1.5 0 0 0 0]}};

test_fun = @(x) x/3.4;                  element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[1/3.4 0 0 0 0]}};
test_fun = @(x) 2.5/x;                  element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[-2.5/(xs_val^2) 0 0 0 0]}};
test_fun = @(x) 5.6./x;                 element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[-5.6/(xs_val^2) 0 0 0 0]}};
test_fun = @(x) x./3.9;                 element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[1/3.9 0 0 0 0]}};
test_fun = @(x) exp(x);                 element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[exp(xs_val) 0 0 0 0]}};

test_fun = @(x) log(x);                 element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[1/xs_val 0 0 0 0]}};
test_fun = @(x) sin(x);                 element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[cos(xs_val) 0 0 0 0]}};
test_fun = @(x) cos(x);                 element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[-sin(xs_val) 0 0 0 0]}};
test_fun = @(x) tan(x);                 element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[sec(xs_val)^2 0 0 0 0]}};
test_fun = @(x) x^3.41;                 element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[3.41*xs_val^2.41 0 0 0 0]}};

test_fun = @(x) x^0;                    element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[0 0 0 0 0]}};
test_fun = @(x) x^(-1);                 element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[-1/(xs_val^2) 0 0 0 0]}};
test_fun = @(x) 3.456^x;                element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[log(3.456)*3.456^xs_val 0 0 0 0]}};
test_fun = @(x) x.^3.41;                element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[3.41*xs_val^2.41 0 0 0 0]}};
test_fun = @(x) x.^0;                   element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[0 0 0 0 0]}};

test_fun = @(x) x.^(-1);                element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[-1/(xs_val^2) 0 0 0 0]}};
test_fun = @(x) 3.456.^x;               element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[log(3.456)*3.456^xs_val 0 0 0 0]}};
test_fun = @(x) sqrt(x);                element1{length(element1)+1}       = {test_fun,{xs},test_fun(xs_val)          ,{[0.5/sqrt(xs_val) 0 0 0 0]}};
%% MATRIX ARITHMETIC OPERATION - SECOND TEST SERIES
%function to test                   Save in as one element                 func      in     corr val                     corr derivative
test_fun = @(x, y, z, w, q) 2.4*x+4.3*y+1.5*z+5.6*w+9.1*q;                 element2{1}                          = {test_fun,{xs ys zs ws qs}, test_fun(xs_val,ys_val,zs_val,ws_val,qs_val), {[2.4 4.3 1.5 5.6 9.1]}};
test_fun = @(x, y, z, w, q) -5.1*x-1.2*y-3.5*z-7.4*w-7.3*q;                element2{length(element2)+1}         = {test_fun,{xs ys zs ws qs}, test_fun(xs_val,ys_val,zs_val,ws_val,qs_val) ,{[-5.1 -1.2 -3.5 -7.4 -7.3]}};
test_fun = @(x, y, z, w, q) 2.4*x-1.2*y-3.5*z+5.6*w-7.3*q;                 element2{length(element2)+1}         = {test_fun,{xs ys zs ws qs}, test_fun(xs_val,ys_val,zs_val,ws_val,qs_val) ,{[2.4 -1.2 -3.5 5.6 -7.3]}};
test_fun = @(x, y, z, w, q) (2.4*x-1.2*y)-3.5*z+(5.6*w-7.3*q);             element2{length(element2)+1}         = {test_fun,{xs ys zs ws qs}, test_fun(xs_val,ys_val,zs_val,ws_val,qs_val) ,{[2.4 -1.2 -3.5 5.6 -7.3]}};
test_fun = @(x, y, z, w, q) (-5.1*x-1.2*y)-(3.5*z-7.4*w-7.3*q);            element2{length(element2)+1}         = {test_fun,{xs ys zs ws qs}, test_fun(xs_val,ys_val,zs_val,ws_val,qs_val) ,{[-5.1 -1.2 -3.5 +7.4 +7.3]}};

test_fun = @(x,y) x*y;                  element2{length(element2)+1}       = {test_fun,{xs ys}, test_fun(xs_val,ys_val) ,{[ys_val xs_val 0 0 0]}};
test_fun = @(x,y) x.*y;                 element2{length(element2)+1}       = {test_fun,{xs ys}, test_fun(xs_val,ys_val) ,{[ys_val xs_val 0 0 0]}};
test_fun = @(x,y) x/y;                  element2{length(element2)+1}       = {test_fun,{xs ys}, test_fun(xs_val,ys_val) ,{[1/ys_val -xs_val/(ys_val^2) 0 0 0]}};
test_fun = @(x,y) x./y;                 element2{length(element2)+1}       = {test_fun,{xs ys}, test_fun(xs_val,ys_val) ,{[1/ys_val -xs_val/(ys_val^2) 0 0 0]}};
test_fun = @(x,y) x./y;                 element2{length(element2)+1}       = {test_fun,{xs ys}, test_fun(xs_val,ys_val) ,{[1/ys_val -xs_val/(ys_val^2) 0 0 0]}};

test_fun = @(x,y) x^y;                  element2{length(element2)+1}       = {test_fun,{xs ys}, test_fun(xs_val,ys_val) ,{[ys_val*xs_val^(ys_val-1) log(xs_val)*xs_val^ys_val 0 0 0]}};
test_fun = @(x,y) x.^y;                 element2{length(element2)+1}       = {test_fun,{xs ys}, test_fun(xs_val,ys_val) ,{[ys_val*xs_val^(ys_val-1) log(xs_val)*xs_val^ys_val 0 0 0]}};
test_fun = @(x,y,z) x + y*z;            element2{length(element2)+1}       = {test_fun,{xs ys zs}, test_fun(xs_val, ys_val, zs_val)         ,{[1 zs_val ys_val 0 0]}};
test_fun = @(x) sin(cos(exp(tan(x))));  element2{length(element2)+1}       = {test_fun,{xs}, test_fun(xs_val)         ,{[-exp(tan(xs_val)) * sec(xs_val)^2 * sin(exp(tan(xs_val))) *cos(cos(exp(tan(xs_val)))) 0 0 0 0 ]}};
%% MATRIX ARITHMETIC OPERATION - THIRD TEST SERIES
test_fun = @(x) x+2.41;                 element3{1}                        = {test_fun,{Bc5}, test_fun(Bc5_val)         ,{eye(5)}};
test_fun = @(x) x-2.81;                 element3{length(element3)+1}       = {test_fun,{Bc5}, test_fun(Bc5_val)         ,{eye(5)}};
test_fun = @(x) 6.41+x;                 element3{length(element3)+1}       = {test_fun,{Bc5}, test_fun(Bc5_val)         ,{eye(5)}};
test_fun = @(x) 6.5 -x;                 element3{length(element3)+1}       = {test_fun,{Bc5}, test_fun(Bc5_val)         ,{-eye(5)}};
test_fun = @(x) -x;                     element3{length(element3)+1}       = {test_fun,{Bc5}, test_fun(Bc5_val)         ,{-eye(5)}};

test_fun = @(x) x*9.2;                  element3{length(element3)+1}       = {test_fun,{Bc5}, test_fun(Bc5_val)         ,{9.2* eye(5)}};
test_fun = @(x) x.*2.11;                element3{length(element3)+1}       = {test_fun,{Bc5}, test_fun(Bc5_val)         ,{2.11*eye(5)}};
test_fun = @(x) 4.2*x;                  element3{length(element3)+1}       = {test_fun,{Bc5}, test_fun(Bc5_val)         ,{4.2* eye(5)}};
test_fun = @(x) 2.46.*x;                element3{length(element3)+1}       = {test_fun,{Bc5}, test_fun(Bc5_val)         ,{2.46*eye(5)}};

test_fun = @(x) x/0.4;                  element3{length(element3)+1}       = {test_fun,{Bc5}, test_fun(Bc5_val)         ,{1/0.4*eye(5)}};
test_fun = @(x) x./2.8;                 element3{length(element3)+1}       = {test_fun,{Bc5}, test_fun(Bc5_val)         ,{1/2.8*eye(5)}};
test_fun = @(x) 5.46./x;                element3{length(element3)+1}       = {test_fun,{Bc5}, test_fun(Bc5_val)         ,{diag(-5.46./(Bc5_val.^2))*eye(5)}};

test_fun = @(x) x.^4.51;                element3{length(element3)+1}       = {test_fun,{Bc5}, test_fun(Bc5_val)         ,{diag(4.51*Bc5_val.^3.51)*eye(5)}};
test_fun = @(x) x.^0;                   element3{length(element3)+1}       = {test_fun,{Bc5}, test_fun(Bc5_val)         ,{0*eye(5)}};
test_fun = @(x) x.^-2.71;               element3{length(element3)+1}       = {test_fun,{Bc5}, test_fun(Bc5_val)         ,{diag(-2.71*Bc5_val.^(-3.71))*eye(5)}};

test_fun = @(x) uminus(x);              element3{length(element3)+1}       = {test_fun,{Bc5}, test_fun(Bc5_val)         ,{-eye(5)}};
test_fun = @(x) exp(x);                 element3{length(element3)+1}       = {test_fun,{Bc5}, test_fun(Bc5_val)         ,{diag(exp(Bc5_val))*eye(5)}};
test_fun = @(x) log(x);                 element3{length(element3)+1}       = {test_fun,{Bc5}, test_fun(Bc5_val)         ,{diag(1./Bc5_val)*eye(5)}};
test_fun = @(x) sin(x);                 element3{length(element3)+1}       = {test_fun,{Bc5}, test_fun(Bc5_val)         ,{diag(cos(Bc5_val))*eye(5)}};
test_fun = @(x) cos(x);                 element3{length(element3)+1}       = {test_fun,{Bc5}, test_fun(Bc5_val)         ,{diag(-sin(Bc5_val))*eye(5)}};
test_fun = @(x) tan(x);                 element3{length(element3)+1}       = {test_fun,{Bc5}, test_fun(Bc5_val)         ,{diag(sec(Bc5_val).^2)*eye(5)}};
test_fun = @(x) sqrt(x);                element3{length(element3)+1}       = {test_fun,{Bc5}, test_fun(Bc5_val)         ,{diag(1/2./sqrt(Bc5_val))*eye(5)}};
%% MATRIX ARITHMETIC OPERATION - FOURTH TEST SERIES
test_fun = @(x) [2 7 3]'.*x;            element4{1}                        = {test_fun,{Cc3}, test_fun(Cc3_val)         ,{[2 0 0;0 7 0;0 0 3]}};
test_fun = @(x) x.*[2 7 3]';            element4{length(element4)+1}       = {test_fun,{Cc3}, test_fun(Cc3_val)         ,{[2 0 0;0 7 0;0 0 3]}};
test_fun = @(x) x./[8;6;2];             element4{length(element4)+1}       = {test_fun,{Cc3}, test_fun(Cc3_val)         ,{[1/8 0 0;0 1/6 0;0 0 1/2]}};
test_fun = @(x) [8 6 2]./x;             element4{length(element4)+1}       = {test_fun,{Cr3}, test_fun(Cr3_val)         ,{-[8 0 0;0 6 0;0 0 2]./diag(Cr3_val.^2)}};

test_fun = @(x) x + [3 7 2]';           element4{length(element4)+1}       = {test_fun,{Cc3}, test_fun(Cc3_val)         ,{eye(3)}};
test_fun = @(x) [3 7 2]' + x;           element4{length(element4)+1}       = {test_fun,{Cc3}, test_fun(Cc3_val)         ,{eye(3)}};

test_fun = @(x) x - [3 7 2]';           element4{length(element4)+1}       = {test_fun,{Cc3}, test_fun(Cc3_val)         ,{eye(3)}};
test_fun = @(x) [3 7 2]' - x;           element4{length(element4)+1}       = {test_fun,{Cc3}, test_fun(Cc3_val)         ,{-eye(3)}};

test_fun = @(x) x.^[3;7;2];             element4{length(element4)+1}       = {test_fun,{Cc3}, test_fun(Cc3_val)         ,{diag([3;7;2].*Cc3.val.^([3;7;2]-1))*eye(3)}};
test_fun = @(x) [3;7;2].^x;             element4{length(element4)+1}       = {test_fun,{Cc3}, test_fun(Cc3_val)         ,{diag([3;7;2].^(Cc3_val).*log([3;7;2]))*eye(3)}};

test_fun = @(x) [4, 3, 7]*x;            element4{length(element4)+1}       = {test_fun,{Cc3}, test_fun(Cc3_val)         ,{[4,3,7]}};
%% MATRIX ARITHMETIC OPERATION - FIFTH TEST SERIES
test_fun = @(x) x.*x;                   element5{1}                        = {test_fun,{Cc3}, test_fun(Cc3_val)         ,{diag(2.*Cc3_val)*eye(3)}};
test_fun = @(x,y) x + y;                element5{length(element5)+1}       = {test_fun,{Cc3 Dc3},test_fun(Cc3_val, Dc3_val) ,{2*eye(3)}};
test_fun = @(x,y) x - y;                element5{length(element5)+1}       = {test_fun,{Cc3 Dc3},test_fun(Cc3_val, Dc3_val) ,{0*eye(3)}};
test_fun = @(x,y) x./ y;                element5{length(element5)+1}       = {test_fun,{Cc3 Dc3},test_fun(Cc3_val, Dc3_val) ,{diag(1./Dc3_val)*eye(3)-diag(Cc3_val./(Dc3_val.^2))}};
%% MATRIX ARITHMETIC OPERATION - SIXTH TEST SERIES PLUS(+)
% AD + VANLIG
test_fun = @(x) x + 21;                 element6{1}                        = {test_fun,{xs} , test_fun(xs_val)          ,{[1 0 0 0 0]              }};
test_fun = @(x) x + [2;7;6];            element6{length(element6)+1}       = {test_fun,{xs} , test_fun(xs_val)          ,{repmat([1 0 0 0 0],3,1)  }};
test_fun = @(x) x + [2 7 6];            element6{length(element6)+1}       = {test_fun,{xs} , test_fun(xs_val)          ,{repmat([1 0 0 0 0]',1,3) }};

test_fun = @(x) x + 21;                 element6{length(element6)+1}       = {test_fun,{Cc3} , test_fun(Cc3_val)       ,{Cc3.der                  }};
test_fun = @(x) x + [2;7;6];            element6{length(element6)+1}       = {test_fun,{Cc3} , test_fun(Cc3_val)       ,{Cc3.der                  }};

test_fun = @(x) x + 23;                 element6{length(element6)+1}       = {test_fun,{Cr3}, test_fun(Cr3_val)        ,{Cr3.der                  }};
test_fun = @(x) x + [2 7 6];            element6{length(element6)+1}       = {test_fun,{Cr3}, test_fun(Cr3_val)        ,{Cr3.der                  }};

% VANLIG + AD
test_fun = @(x) 21 + x;                 element6{length(element6)+1}       = {test_fun,{xs} , test_fun(xs_val)          ,{[1 0 0 0 0]              }};
test_fun = @(x) [2;7;6] + x;            element6{length(element6)+1}       = {test_fun,{xs} , test_fun(xs_val)          ,{repmat([1 0 0 0 0],3,1)  }};
test_fun = @(x) [2 7 6] + x;            element6{length(element6)+1}       = {test_fun,{xs} , test_fun(xs_val)          ,{repmat([1 0 0 0 0]',1,3) }};

test_fun = @(x) 21 + x;                 element6{length(element6)+1}       = {test_fun,{Cc3} , test_fun(Cc3_val)       ,{Cc3.der                  }};
test_fun = @(x) [2;7;6] + x;            element6{length(element6)+1}       = {test_fun,{Cc3} , test_fun(Cc3_val)       ,{Cc3.der                  }};

test_fun = @(x) 23 + x;                 element6{length(element6)+1}       = {test_fun,{Cr3}, test_fun(Cr3_val)        ,{Cr3.der                  }};
test_fun = @(x) [2 7 6] + x;            element6{length(element6)+1}       = {test_fun,{Cr3}, test_fun(Cr3_val)        ,{Cr3.der                  }};

%AD + AD
test_fun = @(x,y) x + y;                element6{length(element6)+1}       = {test_fun,{xs,ys} ,   test_fun(xs_val,ys_val)    ,{[1 1 0 0 0]         }};
test_fun = @(x,y) x + y;                element6{length(element6)+1}       = {test_fun,{xs,Bc5},   test_fun(xs_val,Bc5_val)  ,{eye(5)+[ones(5,1),zeros(5,4)]                }};
test_fun = @(x,y) x + y;                element6{length(element6)+1}       = {test_fun,{xs,Er4},   test_fun(xs_val,Er4_val)  ,{[eye(4);zeros(1,4)] + [ones(1,4);zeros(4,4)] }};

test_fun = @(x,y) x + y;                element6{length(element6)+1}       = {test_fun,{Bc5,xs},   test_fun(xs_val,Bc5_val)  ,{eye(5)+[ones(5,1),zeros(5,4)]                }};
test_fun = @(x,y) x + y;                element6{length(element6)+1}       = {test_fun,{Er4,xs},   test_fun(xs_val,Er4_val)  ,{[eye(4);zeros(1,4)] + [ones(1,4);zeros(4,4)] }};

test_fun = @(x,y) x + y;                element6{length(element6)+1}       = {test_fun,{Bc5,xs},   test_fun(xs_val,Bc5_val)  ,{eye(5)+[ones(5,1),zeros(5,4)]                }};
test_fun = @(x,y) x + y;                element6{length(element6)+1}       = {test_fun,{Er4,xs},   test_fun(xs_val,Er4_val)  ,{[eye(4);zeros(1,4)] + [ones(1,4);zeros(4,4)] }};
%% MATRIX ARITHMETIC OPERATION - SEVENTH TEST SERIES MINUS(-)
% AD + VANLIG
test_fun = @(x) x - 21;                 element7{1}                        = {test_fun,{xs} , test_fun(xs_val)          ,{[1 0 0 0 0]              }};
test_fun = @(x) x - [2;7;6];            element7{length(element7)+1}       = {test_fun,{xs} , test_fun(xs_val)          ,{repmat([1 0 0 0 0],3,1)  }};
test_fun = @(x) x - [2 7 6];            element7{length(element7)+1}       = {test_fun,{xs} , test_fun(xs_val)          ,{repmat([1 0 0 0 0]',1,3) }};

test_fun = @(x) x - 21;                 element7{length(element7)+1}       = {test_fun,{Cc3} , test_fun(Cc3_val)       ,{Cc3.der                  }};
test_fun = @(x) x - [2;7;6];            element7{length(element7)+1}       = {test_fun,{Cc3} , test_fun(Cc3_val)       ,{Cc3.der                  }};

test_fun = @(x) x - 23;                 element7{length(element7)+1}       = {test_fun,{Cr3}, test_fun(Cr3_val)        ,{Cr3.der                  }};
test_fun = @(x) x - [2 7 6];            element7{length(element7)+1}       = {test_fun,{Cr3}, test_fun(Cr3_val)        ,{Cr3.der                  }};

% VANLIG + AD
test_fun = @(x) 21 - x;                 element7{length(element7)+1}       = {test_fun,{xs} , test_fun(xs_val)          ,{-[1 0 0 0 0]              }};
test_fun = @(x) [2;7;6] - x;            element7{length(element7)+1}       = {test_fun,{xs} , test_fun(xs_val)          ,{-repmat([1 0 0 0 0],3,1)  }};
test_fun = @(x) [2 7 6] - x;            element7{length(element7)+1}       = {test_fun,{xs} , test_fun(xs_val)          ,{-repmat([1 0 0 0 0]',1,3) }};

test_fun = @(x) 21 - x;                 element7{length(element7)+1}       = {test_fun,{Cc3} , test_fun(Cc3_val)       ,{-Cc3.der                  }};
test_fun = @(x) [2;7;6] - x;            element7{length(element7)+1}       = {test_fun,{Cc3} , test_fun(Cc3_val)       ,{-Cc3.der                  }};

test_fun = @(x) 23 - x;                 element7{length(element7)+1}       = {test_fun,{Cr3}, test_fun(Cr3_val)        ,{-Cr3.der                  }};
test_fun = @(x) [2 7 6] - x;            element7{length(element7)+1}       = {test_fun,{Cr3}, test_fun(Cr3_val)        ,{-Cr3.der                  }};

%AD + AD
test_fun = @(x,y) x - y;                element7{length(element7)+1}       = {test_fun,{xs,ys} ,   test_fun(xs_val,ys_val)    ,{[1 -1 0 0 0]         }};
test_fun = @(x,y) x - y;                element7{length(element7)+1}       = {test_fun,{xs,Bc5},   test_fun(xs_val,Bc5_val)  ,{-eye(5)+[ones(5,1),zeros(5,4)]                }};
test_fun = @(x,y) x - y;                element7{length(element7)+1}       = {test_fun,{xs,Er4},   test_fun(xs_val,Er4_val)  ,{-[eye(4);zeros(1,4)] + [ones(1,4);zeros(4,4)] }};

test_fun = @(x,y) x - y;                element7{length(element7)+1}       = {test_fun,{Bc5,xs},   test_fun(Bc5_val,xs_val)  ,{eye(5)-[ones(5,1),zeros(5,4)]                 }};
test_fun = @(x,y) x - y;                element7{length(element7)+1}       = {test_fun,{Bc5,Fc5},  test_fun(Bc5_val,Fc5_val) ,{zeros(5)                                      }};

test_fun = @(x,y) x - y;                element7{length(element7)+1}       = {test_fun,{Er4,xs} ,  test_fun(Er4_val,xs_val)  ,{[eye(4);zeros(1,4)]-[ones(1,4);zeros(4,4)]    }};
test_fun = @(x,y) x - y;                element7{length(element7)+1}       = {test_fun,{Er4,Gr4},  test_fun(Er4_val,Gr4_val) ,{zeros(5,4)                                    }};
%% MATRIX ARITHMETIC OPERATION - EIGHT TEST SERIES TIMES(.*)
% AD + VANLIG
test_fun = @(x) x .* 31;                 element8{1}                       = {test_fun,{xs}, test_fun(xs_val)          ,{[31 0 0 0 0 ]                           }};
test_fun = @(x) x .* [2;7;6];            element8{length(element8)+1}      = {test_fun,{xs}, test_fun(xs_val)          ,{[[2 0 0 0 0];[7 0 0 0 0];[6 0 0 0 0]]   }};
test_fun = @(x) x .* [2 7 6];            element8{length(element8)+1}      = {test_fun,{xs}, test_fun(xs_val)          ,{[[2 0 0 0 0]',[7 0 0 0 0]',[6 0 0 0 0]']}};

test_fun = @(x) x .* 8;                  element8{length(element8)+1}      = {test_fun,{Bc5}, test_fun(Bc5_val)        ,{8*eye(5)                                }};
test_fun = @(x) x .* [2;7;6;8;6];        element8{length(element8)+1}      = {test_fun,{Bc5}, test_fun(Bc5_val)        ,{diag([2;7;6;8;6])                       }};

test_fun = @(x) x .* 5;                  element8{length(element8)+1}      = {test_fun,{Er4}, test_fun(Er4_val)        ,{5*[eye(4);zeros(1,4)]                   }};
test_fun = @(x) x .* [1 4 2 3];          element8{length(element8)+1}      = {test_fun,{Er4}, test_fun(Er4_val)        ,{  [diag([1 4 2 3]);zeros(1,4)]          }};

% VANLIG + AD
test_fun = @(x) 31 .* x;                 element8{length(element8)+1}      = {test_fun,{xs}, test_fun(xs_val)          ,{[31 0 0 0 0 ]                           }};
test_fun = @(x) [2;7;6] .* x;            element8{length(element8)+1}      = {test_fun,{xs}, test_fun(xs_val)          ,{[[2 0 0 0 0];[7 0 0 0 0];[6 0 0 0 0]]   }};
test_fun = @(x) [2 7 6] .* x;            element8{length(element8)+1}      = {test_fun,{xs}, test_fun(xs_val)          ,{[[2 0 0 0 0]',[7 0 0 0 0]',[6 0 0 0 0]']}};

test_fun = @(x) 8 .* x;                  element8{length(element8)+1}      = {test_fun,{Bc5}, test_fun(Bc5_val)        ,{8*eye(5)                                }};
test_fun = @(x) [2;7;6;8;6] .* x;        element8{length(element8)+1}      = {test_fun,{Bc5}, test_fun(Bc5_val)        ,{diag([2;7;6;8;6])                       }};

test_fun = @(x) 5 .* x;                  element8{length(element8)+1}      = {test_fun,{Er4}, test_fun(Er4_val)        ,{5*[eye(4);zeros(1,4)]                   }};
test_fun = @(x) [1 4 2 3] .* x;          element8{length(element8)+1}      = {test_fun,{Er4}, test_fun(Er4_val)        ,{  [diag([1 4 2 3]);zeros(1,4)]          }};

%AD + AD
test_fun = @(x,y) x .* y;                element8{length(element8)+1}      = {test_fun,{xs,ys} ,   test_fun(xs_val,ys_val)   ,{[ys_val xs_val 0 0 0]                  }};
test_fun = @(x,y) x .* y;                element8{length(element8)+1}      = {test_fun,{xs,Bc5},   test_fun(xs_val,Bc5_val)  ,{xs_val*eye(5) + [ones(5,1).*Bc5_val,zeros(5,4)] }};
test_fun = @(x,y) x .* y;                element8{length(element8)+1}      = {test_fun,{xs,Er4},   test_fun(xs_val,Er4_val)  ,{xs_val*[eye(4);zeros(1,4)] + [Er4_val.*ones(1,4);zeros(4)] }};

test_fun = @(x,y) x .* y;                element8{length(element8)+1}      = {test_fun,{Bc5,xs},   test_fun(Bc5_val,xs_val)  ,{xs_val*eye(5) + [ones(5,1).*Bc5_val,zeros(5,4)] }};
test_fun = @(x,y) x .* y;                element8{length(element8)+1}      = {test_fun,{Bc5,Fc5},  test_fun(Bc5_val,Fc5_val) ,{diag(Bc5_val) + diag(Fc5_val)       }};

test_fun = @(x,y) x .* y;                element8{length(element8)+1}      = {test_fun,{Er4,xs},   test_fun(xs_val,Er4_val)  ,{xs_val*[eye(4);zeros(1,4)] + [ones(1,4).*Er4_val;zeros(4)]  }};
test_fun = @(x,y) x .* y;                element8{length(element8)+1}      = {test_fun,{Er4,Gr4},  test_fun(Er4_val,Gr4_val) ,{[[diag(Er4_val) + diag(Gr4_val)];zeros(1,4)]            }};
%% MATRIX ARITHMETIC OPERATION - NINTH TEST SERIES POWER(.^)
% AD + VANLIG
test_fun = @(x) x .^ 6;                  element9{1}                       = {test_fun,{xs}, test_fun(xs_val)          ,{[6.*xs_val.^5 0 0 0 0]                           }};
test_fun = @(x) x .^ [2;7;6];            element9{length(element9)+1}      = {test_fun,{xs}, test_fun(xs_val)          ,{[([2;7;6]).*xs_val.^[1;6;5],zeros(3,4)]          }};
test_fun = @(x) x .^ [2 7 6];            element9{length(element9)+1}      = {test_fun,{xs}, test_fun(xs_val)          ,{[([2 7 6]).*xs_val.^[1 6 5];zeros(4,3)]          }};

test_fun = @(x) x .^ 6;                  element9{length(element9)+1}      = {test_fun,{Bc5}, test_fun(Bc5_val)        ,{diag(6.*Bc5_val.^5)*eye(5)                       }};
test_fun = @(x) x .^ [2;7;6;1;4];        element9{length(element9)+1}      = {test_fun,{Bc5}, test_fun(Bc5_val)        ,{diag([2;7;6;1;4].*Bc5_val.^[1;6;5;0;3])          }};

test_fun = @(x) x .^ 6;                  element9{length(element9)+1}      = {test_fun,{Er4}, test_fun(Er4_val)        ,{[diag(6.*Er4_val.^5)*eye(4);zeros(1,4)]          }};
test_fun = @(x) x .^ [5 3 2 3];          element9{length(element9)+1}      = {test_fun,{Er4}, test_fun(Er4_val)        ,{[diag([5 3 2 3].*Er4_val.^[4 2 1 2]);zeros(1,4)] }};

% VANLIG + AD 
test_fun = @(x) 6 .^ x;                  element9{length(element9)+1}      = {test_fun,{xs}, test_fun(xs_val)          ,{[6.^xs_val.*log(6) 0 0 0 0]                           }};
test_fun = @(x) [2;7;6] .^ x;            element9{length(element9)+1}      = {test_fun,{xs}, test_fun(xs_val)          ,{[[2;7;6].^xs_val.*log([2;7;6]).*ones(3,1),zeros(3,4)]  }};
test_fun = @(x) [2 7 6] .^ x;            element9{length(element9)+1}      = {test_fun,{xs}, test_fun(xs_val)          ,{[[2;7;6].^xs_val.*log([2;7;6]).*ones(3,1),zeros(3,4)]'  }};

test_fun = @(x) 6 .^ x;                  element9{length(element9)+1}      = {test_fun,{Bc5}, test_fun(Bc5_val)        ,{6.^Bc5_val.*log(6).*eye(5)          }};
test_fun = @(x) [2;7;6;1;4] .^ x;        element9{length(element9)+1}      = {test_fun,{Bc5}, test_fun(Bc5_val)        ,{diag([2;7;6;1;4].^Bc5_val.*log([2;7;6;1;4]))*eye(5)          }};

test_fun = @(x) 6 .^ x;                  element9{length(element9)+1}      = {test_fun,{Er4}, test_fun(Er4_val)        ,{[diag(6.^Er4_val.*log(6));zeros(1,4)]         }};
test_fun = @(x) [5 3 2 3] .^ x;          element9{length(element9)+1}      = {test_fun,{Er4}, test_fun(Er4_val)        ,{ [diag([5 3 2 3].^Er4_val.*log([5 3 2 3]));zeros(1,4)] }};
%% MATRIX ARITHMETIC OPERATION - TENTH TEST SERIES RDIVIDE(./)
test_fun = @(x) x ./ 4;                  element10{1}                      = {test_fun,{xs}, test_fun(xs_val)          ,{[1/4 0 0 0 0]                             }};
test_fun = @(x) x ./ [2;7;6];            element10{length(element10)+1}    = {test_fun,{xs}, test_fun(xs_val)          ,{[[1/2;1/7;1/6],zeros(3,4)]                }};
test_fun = @(x) x ./ [2 7 6];            element10{length(element10)+1}    = {test_fun,{xs}, test_fun(xs_val)          ,{[[1/2 1/7 1/6];zeros(4,3)]                }};

test_fun = @(x) x ./ 4;                  element10{length(element10)+1}    = {test_fun,{Bc5}, test_fun(Bc5_val)        ,{1/4*eye(5)                                }};
test_fun = @(x) x ./ [2;7;6;3;8];        element10{length(element10)+1}    = {test_fun,{Bc5}, test_fun(Bc5_val)        ,{diag(1./[2;7;6;3;8])                      }};

test_fun = @(x) x ./ 7;                  element10{length(element10)+1}    = {test_fun,{Er4}, test_fun(Er4_val)        ,{1/7*[eye(4);zeros(1,4)]                   }};
test_fun = @(x) x ./ [2 7 6 1];          element10{length(element10)+1}    = {test_fun,{Er4}, test_fun(Er4_val)        ,{ [diag(1./[2 7 6 1]);zeros(1,4)]          }};

test_fun = @(x) 4 ./ x;                  element10{length(element10)+1}    = {test_fun,{xs}, test_fun(xs_val)          ,{[-4/xs_val^2 0 0 0 0]                     }};
test_fun = @(x) [2;7;6] ./ x;            element10{length(element10)+1}    = {test_fun,{xs}, test_fun(xs_val)          ,{[-[2;7;6]./xs_val^2,zeros(3,4)]           }};
test_fun = @(x) [2 7 6] ./ x;            element10{length(element10)+1}    = {test_fun,{xs}, test_fun(xs_val)          ,{[-[2;7;6]./xs_val^2,zeros(3,4)]'          }};

test_fun = @(x) 4 ./ x;                  element10{length(element10)+1}    = {test_fun,{Bc5}, test_fun(Bc5_val)        ,{ diag(-4./(Bc5_val.^2))                   }};
test_fun = @(x) [2;7;6;3;8] ./ x;        element10{length(element10)+1}    = {test_fun,{Bc5}, test_fun(Bc5_val)        ,{ diag(-[2;7;6;3;8]./(Bc5_val.^2))         }};

test_fun = @(x) 7 ./ x;                  element10{length(element10)+1}    = {test_fun,{Er4}, test_fun(Er4_val)        ,{[diag(-7./Er4_val.^2);zeros(1,4)]         }};
test_fun = @(x) [2 7 6 1] ./ x;          element10{length(element10)+1}    = {test_fun,{Er4}, test_fun(Er4_val)        ,{[diag(-[2 7 6 1]./Er4_val.^2);zeros(1,4)] }};

%AD + AD
test_fun = @(x,y) x ./ y;                element10{length(element10)+1}      = {test_fun,{xs,ys} ,   test_fun(xs_val,ys_val)   ,{[1/ys_val -xs_val./ys_val^2 0 0 0]                  }};
test_fun = @(x,y) x ./ y;                element10{length(element10)+1}      = {test_fun,{xs,Bc5},   test_fun(xs_val,Bc5_val)  ,{[1./Bc5_val,zeros(5,4)] - diag(xs_val./Bc5_val.^2) }};
test_fun = @(x,y) x ./ y;                element10{length(element10)+1}      = {test_fun,{xs,Er4},   test_fun(xs_val,Er4_val)  ,{[1./Er4_val;zeros(4,4)] - [diag(xs_val./Er4_val.^2);zeros(1,4)] }};

test_fun = @(x,y) x ./ y;                element10{length(element10)+1}      = {test_fun,{Bc5,xs},   test_fun(Bc5_val,xs_val)  ,{-[Bc5_val./xs_val^2,zeros(5,4)] + eye(5)./xs_val }};
test_fun = @(x,y) x ./ y;                element10{length(element10)+1}      = {test_fun,{Bc5,Fc5},  test_fun(Bc5_val,Fc5_val) ,{-diag(Bc5_val./Fc5_val.^2) + eye(5)./Fc5_val }};

test_fun = @(x,y) x ./ y;                element10{length(element10)+1}      = {test_fun,{Er4,xs},   test_fun(Er4_val,xs_val)  ,{-[Er4_val./xs_val^2;zeros(4,4)] + [eye(4);zeros(1,4)]./xs_val  }};
test_fun = @(x,y) x ./ y;                element10{length(element10)+1}      = {test_fun,{Er4,Gr4},  test_fun(Er4_val,Gr4_val) ,{-[diag(Er4_val./Gr4_val.^2);zeros(1,4)] + [diag(1./Gr4_val);zeros(1,4)]            }};
%% MATRIX ARITHMETIC OPERATION - ELEVENTH TEST SERIES MRDIVIDE(/)
% AD + VANLIG
test_fun = @(x) x / 4;                  element11{1}                      = {test_fun,{xs}, test_fun(xs_val)          ,{[1/4 0 0 0 0]                             }};
test_fun = @(x) x / 4;                  element11{length(element11)+1}    = {test_fun,{Bc5}, test_fun(Bc5_val)        ,{1/4*eye(5)                                }};
test_fun = @(x) x / 7;                  element11{length(element11)+1}    = {test_fun,{Er4}, test_fun(Er4_val)        ,{1/7*[eye(4);zeros(1,4)]                   }};

%VANLIG + AD
test_fun = @(x) 4 / x;                  element11{length(element11)+1}    = {test_fun,{xs}, test_fun(xs_val)          ,{[-4/xs_val^2 0 0 0 0]                     }};
test_fun = @(x) [2;7;6] / x;            element11{length(element11)+1}    = {test_fun,{xs}, test_fun(xs_val)          ,{[-[2;7;6]./xs_val^2,zeros(3,4)]           }};
test_fun = @(x) [2 7 6] / x;            element11{length(element11)+1}    = {test_fun,{xs}, test_fun(xs_val)          ,{[-[2;7;6]./xs_val^2,zeros(3,4)]'          }};

%AD + AD
test_fun = @(x,y) x / y;                element11{length(element11)+1}    = {test_fun,{xs,ys} ,   test_fun(xs_val,ys_val)   ,{[1/ys_val -xs_val./ys_val^2 0 0 0]                  }};
test_fun = @(x,y) x / y;                element11{length(element11)+1}    = {test_fun,{Bc5,xs},   test_fun(Bc5_val,xs_val)  ,{-[Bc5_val./xs_val^2,zeros(5,4)] + eye(5)./xs_val }};
test_fun = @(x,y) x / y;                element11{length(element11)+1}    = {test_fun,{Er4,xs},   test_fun(Er4_val,xs_val)  ,{-[Er4_val./xs_val^2;zeros(4,4)] + [eye(4);zeros(1,4)]./xs_val  }};
%% MATRIX ARITHMETIC OPERATION - TWELTH TEST SERIES MPOWER(^)
test_fun = @(x) x ^ 6;                   element12{1}                      = {test_fun,{xs}, test_fun(xs_val)             ,{[6.*xs_val.^5 0 0 0 0]}};
test_fun = @(x) 6 ^ x;                   element12{length(element12)+1}    = {test_fun,{xs}, test_fun(xs_val)             ,{[log(6).*6.^xs_val 0 0 0 0]}};
test_fun = @(x,y) x ^ y;                 element12{length(element12)+1}    = {test_fun,{xs,ys},   test_fun(xs_val,ys_val)   ,{[ys_val.*xs_val.^(ys_val-1) log(xs_val).*xs_val.^ys_val 0 0 0]}};
%% MATRIX ARITHMETIC OPERATION - THIRTEENTH TEST SERIES MTIMES(*)
% AD + VANLIG
test_fun = @(x) x * 31;                 element13{1}                       = {test_fun,{xs}, test_fun(xs_val)          ,{[31 0 0 0 0 ]                           }};
test_fun = @(x) x * [2;7;6];            element13{length(element13)+1}      = {test_fun,{xs}, test_fun(xs_val)          ,{[[2 0 0 0 0];[7 0 0 0 0];[6 0 0 0 0]]   }};
test_fun = @(x) x * [2 7 6];            element13{length(element13)+1}      = {test_fun,{xs}, test_fun(xs_val)          ,{[[2 0 0 0 0]',[7 0 0 0 0]',[6 0 0 0 0]']}};

test_fun = @(x) x * 8;                  element13{length(element13)+1}      = {test_fun,{Bc5}, test_fun(Bc5_val)        ,{8*eye(5)                                }};
test_fun = @(x) x * 5;                  element13{length(element13)+1}      = {test_fun,{Er4}, test_fun(Er4_val)        ,{5*[eye(4);zeros(1,4)]                   }};

% VANLIG + AD
test_fun = @(x) 31 * x;                 element13{length(element13)+1}      = {test_fun,{xs}, test_fun(xs_val)          ,{[31 0 0 0 0 ]                           }};
test_fun = @(x) [2;7;6] * x;            element13{length(element13)+1}      = {test_fun,{xs}, test_fun(xs_val)          ,{[[2 0 0 0 0];[7 0 0 0 0];[6 0 0 0 0]]   }};
test_fun = @(x) [2 7 6] * x;            element13{length(element13)+1}      = {test_fun,{xs}, test_fun(xs_val)          ,{[[2 0 0 0 0]',[7 0 0 0 0]',[6 0 0 0 0]']}};

test_fun = @(x) 8 * x;                  element13{length(element13)+1}      = {test_fun,{Bc5}, test_fun(Bc5_val)        ,{8*eye(5)                                }};
test_fun = @(x) 5 * x;                  element13{length(element13)+1}      = {test_fun,{Er4}, test_fun(Er4_val)        ,{5*[eye(4);zeros(1,4)]                   }};

%AD + AD
test_fun = @(x,y) x * y;                element13{length(element13)+1}      = {test_fun,{xs,ys} ,   test_fun(xs_val,ys_val)   ,{[ys_val xs_val 0 0 0]                  }};
test_fun = @(x,y) x * y;                element13{length(element13)+1}      = {test_fun,{xs,Bc5},   test_fun(xs_val,Bc5_val)  ,{xs_val*eye(5) + [ones(5,1).*Bc5_val,zeros(5,4)] }};
test_fun = @(x,y) x * y;                element13{length(element13)+1}      = {test_fun,{xs,Er4},   test_fun(xs_val,Er4_val)  ,{xs_val*[eye(4);zeros(1,4)] + [Er4_val.*ones(1,4);zeros(4)] }};

test_fun = @(x,y) x * y;                element13{length(element13)+1}      = {test_fun,{Bc5,xs},   test_fun(Bc5_val,xs_val)  ,{xs_val*eye(5) + [ones(5,1).*Bc5_val,zeros(5,4)] }};
test_fun = @(x,y) x * y;                element13{length(element13)+1}      = {test_fun,{Er4,xs},   test_fun(xs_val,Er4_val)  ,{xs_val*[eye(4);zeros(1,4)] + [ones(1,4).*Er4_val;zeros(4)]  }};

%% MATRIX ARITHMETIC OPERATION - TWELTH TEST SERIES MPOWER(^)
test_fun = @(x) normcdf(x,0,1);         element14{1}                        = {test_fun,{xs}, test_fun(xs_val)             ,{[normpdf(xs_val,0,1) 0 0 0 0]}};

%% MATRIX SHAPE      OPERATION - FIRST TEST SERIE
% first = 8;          second = 9;         op_element{1}                      = {'two',first,second,I_val(first,second),I_id(first,second)};
% first = 12;         second = 14;        op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
% first = 2;          second = 4;         op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
% first = ':';        second = 4;         op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
% first = 6;          second = ':';       op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
%
% first = ':';        second = ':';       op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
% first = 2;          second = 2:5;       op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
% first = 4;          second = 2:2:15;    op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
% first = 3:7;        second = 6;         op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
% first = 3:2:9;      second = 9;         op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
%
% first = ':';        second = 2:8;       op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
% first = ':';        second = 3:4:18;    op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
% first = 3:9;        second = ':';       op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
% first = 4:6:20;     second = ':';       op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
% first = 6:15;       second = 8:12;      op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
%
% first = 4:6:20;     second = 3:4:18;    op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
% first = 4;          second = [];        op_element{length(op_element)+1}   = {'one',first,second,I_val(first),I_id(first)};
% first = 55;         second = [];        op_element{length(op_element)+1}   = {'one',first,second,I_val(first),I_id(first)};
% first = 45:55;      second = [];        op_element{length(op_element)+1}   = {'one',first,second,I_val(first),I_id(first)};
% first = 3:5:68;     second = [];        op_element{length(op_element)+1}   = {'one',first,second,I_val(first),I_id(first)};
%
% first = ':';        second = [];        op_element{length(op_element)+1}   = {'one',first,second,I_val(first),I_id(first)};
%% ===========OUTPUT=======================================================
%% ARITHMETIC OPERATION - VECTOR TO TEST
n_1 = length(element1);
n_2 = length(element2);
n_3 = length(element3);
n_4 = length(element4);
n_5 = length(element5);
n_6 = length(element6);
n_7 = length(element7);
n_8 = length(element8);
n_9 = length(element9);
n_10 = length(element10);
n_11 = length(element11);
n_12 = length(element12);
n_13 = length(element13);
n_14 = length(element14);

name_1 = {'FIRST - SCALAR ONE DIMENSIONAL'};
name_2 = {'SECOND - SCALAR TWO DIMENSIONAL'};
name_3 = {'THIRD - MATRIX ONE DIMENSIONAL'};
name_4 = {'FOURTH - MATRIX TWO DIMENSIONAL'};
name_5 = {'FIFTH - HARD COMBINATION'};
name_6 = {'SIXTH - REST COMBINATION'};
name_7 = {'SEVENTH - REST COMBINATION'};
name_8 = {'EIGHT - REST COMBINATION'};
name_9 = {'NINTH - REST COMBINATION'};
name_10= {'TENTH - REST COMBINATION'};
name_11= {'ELEVENTH - REST COMBINATION'};
name_12= {'TWELTH - REST COMBINATION'};
name_13= {'THIRTEENTH - REST COMBINATION'};
name_14= {'FOURTEENTH - PROBABILITY'};

n = [n_1 n_2 n_3 n_4 n_5 n_6 n_7 n_8 n_9 n_10 n_11 n_12 n_13 n_14];
name = {name_1, name_2, name_3,name_4,name_5,name_6,name_7,name_8,name_9,name_10,name_11,name_12,name_13,name_14};
element = {element1,element2,element3,element4,element5,element6,element7,element8,element9,element10,element11,element12,element13,element14};
%
test_vector = get_tv(n,name,element);
%% ARITHMETIC OPERATION TEST ALGORITHM
n_t_vec = length(test_vector);
i_test = 0;
tic
% LOOP ALL TESTS IN THE SERIES
for i=1:n_t_vec
    %READ IN
    temp_str            = test_vector{i};
    try
        if(~ischar(temp_str))
            %PRE-PROCESSING
            i_test = i_test + 1;
            temp_fun            = temp_str{1};
            temp_eval_pnt_str   = temp_str{2};
            temp_corr_val       = temp_str{3};
            temp_corr_der_str   = temp_str{4};
%             temp_n_in           = temp_str{5};
            
            %% EVALUATION OF FUNCTION VALUE AND DERIVATIVE
            if(iscell(temp_eval_pnt_str))
                if(5==length(temp_eval_pnt_str))
                    temp_eval = feval(temp_fun,temp_eval_pnt_str{1,1},temp_eval_pnt_str{1,2},temp_eval_pnt_str{1,3},temp_eval_pnt_str{1,4},temp_eval_pnt_str{1,5});
                elseif(4==length(temp_eval_pnt_str))
                    temp_eval = feval(temp_fun,temp_eval_pnt_str{1,1},temp_eval_pnt_str{1,2},temp_eval_pnt_str{1,3},temp_eval_pnt_str{1,4});
                elseif(3==length(temp_eval_pnt_str))
                    temp_eval = feval(temp_fun,temp_eval_pnt_str{1,1},temp_eval_pnt_str{1,2},temp_eval_pnt_str{1,3});
                elseif(2==length(temp_eval_pnt_str))
                    temp_eval = feval(temp_fun,temp_eval_pnt_str{1,1},temp_eval_pnt_str{1,2});
                elseif(1==length(temp_eval_pnt_str))
                    temp_eval = feval(temp_fun,temp_eval_pnt_str{1,1});
                else
                    error('error');
                end
            else
                temp_eval = feval(temp_fun,temp_eval_pnt_str);
            end
            [n_row,n_col] = size(temp_eval);
            val_bool = true;
            der_bool = true;
            
            Bt = temp_eval.der;
            %Value
            if(abs(temp_corr_val-temp_eval.val)<10^(-14))
                val_bool = val_bool*true;
            else
                val_bool = val_bool*false;
            end
            temp_corr_der = temp_corr_der_str{1,1};
            %Derivative
            if(max(abs(Bt - temp_corr_der))<10^(-14))
                der_bool = der_bool*true;
            else
                der_bool = der_bool*false;
            end
            
            %% OUTPUT
            %Dela upp outputen i delar av fem
            if(mod(i_test-1,5)==0 && i_test~=1)
                fprintf(fid1,'\n');
            end
            if(val_bool)
                if(i_test>=10)
                    fprintf(fid1,'The value for test %.0f is \t    correct ',i_test);
                else
                    fprintf(fid1,'The value for test %.0f  is \t    correct ',i_test);
                end
            else
                fprintf(fid1,'The value for test %0.f is \t NOT correct ',i_test);
            end
            if(der_bool)
                if(i_test>=10)
                    fprintf(fid1,'and the derivatives for test %0.f is \t           correct.\n',i_test);
                else
                    fprintf(fid1,'and the derivatives for test %0.f  is \t           correct.\n',i_test);
                end
            else
                if(i_test>=10)
                    fprintf(fid1,'and the derivatives for test %0.f is \tNOT        correct        %-4.0f\n',i_test,i);
                else
                    fprintf(fid1,'and the derivatives for test %0.f  is \tNOT        correct       %-4.0f\n',i_test,i);
                end
            end
        else
            fprintf(fid1,'===================%-30s=============================================\n',temp_str);
            i_test = 0;
        end
    catch me
        fprintf(fid1,    '*******************ERROR FOR TEST %-4.0f***************************************************************%-4.0f',i_test,i);
        if(strcmp(me.identifier,'MATLAB:sub2ind:IndexOutOfRange'))
            fprintf(fid1,'CHECK The NUMBER OF INPUT VARIABLES'); %This is a known error in the unit test but not in the code
        else
            fprintf(fid1,me.identifier);
        end
        fprintf(fid1,'\n');
    end
end
toc

%AFTER PROCESSING FOR THE FIRST SERIES
fclose(fid1);
%% SHAPE      OPERATION TEST ALGORITHM

% for i_test=1:length(op_element)
%
%     t_element = op_element{i_test};
%     op_case   = t_element{1};
%     first     = t_element{2};
%     second    = t_element{3};
%     r_val     = t_element{4};
%     r_id      = t_element{5};
%
%     switch op_case
%         case 'one'
%             b = I(first);
%             if(b.val==r_val)
%                 bool_val = true;
%             else
%                 bool_val = false;
%             end
%             if(b.id==r_id)
%                 bool_id = true;
%             else
%                 bool_id = false;
%             end
%         case 'two'
%             b = I(first,second);
%             if(b.val==r_val)
%                 bool_val = true;
%             else
%                 bool_val = false;
%             end
%             if(b.id==r_id)
%                 bool_id = true;
%             else
%                 bool_id = false;
%             end
%     end
%
%     if(mod(i_test-1,5)==0 && i_test~=1)
%         fprintf(fid2,'\n');
%     end
%
%     if(bool_val)
%         if(i_test>=10)
%             fprintf(fid2,'The value for test %0.f is \t    correct ',i_test);
%         else
%             fprintf(fid2,'The value for test %0.f  is \t    correct ',i_test);
%         end
%     else
%         fprintf(fid2,'The value for test %0.f is \t NOT correct ',i_test);
%     end
%
%     if(bool_id)
%         if(i_test>=10)
%             fprintf(fid2,'and the id for test %0.f is \t           correct.\n',i_test);
%         else
%             fprintf(fid2,'and the id for test %0.f  is \t           correct.\n',i_test);
%         end
%     else
%         if(i_test>=10)
%             fprintf(fid2,'and the id for test %0.f is \tNOT        correct.\n',i_test);
%         else
%             fprintf(fid2,'and the id for test %0.f  is \tNOT        correct.\n',i_test);
%         end
%     end
%
% end

fclose(fid2);