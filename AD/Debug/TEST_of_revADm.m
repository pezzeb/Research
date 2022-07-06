
%% INIT-LOADING

clear;
clc;
clear global tape;
global tape;
global tape_switch;

addpath('.\class');
addpath('..\function');
addpath('.\test_function');

init_size_tape = 10000; 
tape_switch = 'main_tape';

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
fid1 = fopen('.\outputRev\output_from_TEST_of_revADm.txt','w');
fid2 = fopen('.\outputRev\output_from_TEST_of_revADm_matrix_operations.txt','w');
%**************************************************************************

%% INPUT VARIABLES - OPERATIONS

n_in = 9;
tape = -speye(n_in);

x_val = 3;
y_val = 2;
z_val = 3;
w_val = 6;
q_val = 12;
j_val = 8;

reset_tape(5);

x = revADm(x_val,1);
y = revADm(y_val,2);
z = revADm(z_val,3);
w = revADm(w_val,4);
q = revADm(q_val,5);
j = revADm(j_val,0);

yy_val = [8,9];
yy_id  = [2 3];
yy = revADm(yy_val,yy_id);

%Constants
A_val = [3 2;4 5];
A_id  = [1 3;2 4];
n_in_A = 4;
a1_val = 3;
a2_val = 4;
a3_val = 2;
a4_val = 5;

AA1_val = [8 3;4 5];
AA1_id  = [1 0;0 2];

AA2_val = [2 1;6 2];
AA2_id  = [0 3;0 0];

B_val = [1.1 5 6;2 3 4];
Bs_val= [0.1 0.2 0.3;0.4 0.5 0.6];
B_id  = [1 3 5;2 4 6];
n_in_B= 6;

C_val = 0.2*[1 5 6;2 3 4];
C_id  = [7 9 11;8 10 12];
n_in_C= 6;

E_val = [5 6;2 3;1 5;6 9];
E_id  = [1 5;2 6;3 7;4 8];
n_in_E = 8;

F_val   = rand(10,10);
F_id    = reshape(1:1:100,10,10);
n_in_F  = 10^2;
Fb      = rand(10,1);
FAb     = rand(15,10);
FAa     = rand(10,18);

% row_G = 20;
% col_G = 120;
% G_val   = rand(row_G,col_G);
% G_id    = reshape(1:1:numel(G_val),row_G,col_G);
% n_in_G  = numel(G_val);
% 
% row_H = col_G;
% col_H = 42;
% H_val   = rand(row_H,col_H);
% H_id    = reshape(1:1:numel(H_val),row_H,col_H)+n_in_G;
% n_in_H  = numel(H_val);

row_G = 20;
col_G = 120;
G_val   = rand(row_G,col_G);
G_id    = reshape(1:1:numel(G_val),row_G,col_G);
n_in_G  = numel(G_val);

row_H = col_G;
col_H = 42;
H_val   = rand(row_H,col_H);
H_id    = reshape(1:1:numel(H_val),row_H,col_H)+n_in_G;
n_in_H  = numel(H_val);


b_val = [4 5 6 7 8 9]; bs_val = [0.1 0.2 0.3 0.4 0.5 0.6];
b_id  = [1 2 3 4 5 6];
n_in_b = 6;

c_val = [4 5 6 7];
c_id  = [1 2 3 4];
n_in_c = 4;

d_val = [3 1 4 9];
d_id  = [5 6 7 8];
n_in_d = 4;

e_val = [6 2];
e_id  = [5 6];
n_in_e = 2;

%AD variables
A = revADm(A_val,A_id);
AA1= revADm(AA1_val,AA1_id);
AA2= revADm(AA2_val,AA2_id);
b = revADm(b_val,b_id);bs = revADm(bs_val,b_id); %Den andra används om det finns begränsningar
B = revADm(B_val,B_id);Bs = revADm(Bs_val,B_id); %Den andra används om det finns begränsningar
C = revADm(C_val,C_id);
c = revADm(c_val,c_id);
d = revADm(d_val,d_id);
E = revADm(E_val,E_id);
e = revADm(e_val,e_id);
F = revADm(F_val,F_id);
G = revADm(G_val,G_id);
H = revADm(H_val,H_id);
%I används nedan


reset_tape(n_in_F);
% % C = [3 5]*B;
% reset_tape(8);
% D = F*Fb;

% l = A(:);
reset_tape(5);
f = 2.4*x+5;

%% INPUT VARIABLES - MATRIX OPERATIONS
    n_r_I = 25;
n_c_I = 65;

I_val = rand(n_r_I,n_c_I);
I_id  = reshape(1:1:numel(I_val),n_r_I,n_c_I);
I     = revADm(I_val,I_id);

%% MATRIX ARITHMETIC OPERATION - FIRST TEST SERIES
%%SCALAR TEST
%function to test                   Save in as one element                 func     in  corr val          corr derivative
test_fun = @(x) x + 2;                  element1{1}                        = {test_fun,{x},test_fun(x_val)          ,{[1 0 0 0 0]},5};
test_fun = @(x) x + 6.42;               element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[1 0 0 0 0]},5};
test_fun = @(x) x - 6.42;               element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[1 0 0 0 0]},5};
test_fun = @(x) 6.9 - x;                element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[-1 0 0 0 0]},5};
test_fun = @(x) 6.9 + x;                element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[1 0 0 0 0]},5};

test_fun = @(x) -x;                     element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[-1 0 0 0 0]},5};
test_fun = @(x) 5.5*x;                  element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[5.5 0 0 0 0]},5};
test_fun = @(x) x*3.4;                  element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[3.4 0 0 0 0]},5};
test_fun = @(x) 7.2.*x;                 element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[7.2 0 0 0 0]},5};
test_fun = @(x) x.*1.5;                 element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[1.5 0 0 0 0]},5};

test_fun = @(x) x/3.4;                  element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[1/3.4 0 0 0 0]},5};
test_fun = @(x) 2.5/x;                  element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[-2.5/(x_val^2) 0 0 0 0]},5};
test_fun = @(x) 5.6./x;                 element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[-5.6/(x_val^2) 0 0 0 0]},5};
test_fun = @(x) x./3.9;                 element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[1/3.9 0 0 0 0]},5};
test_fun = @(x) exp(x);                 element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[exp(x_val) 0 0 0 0]},5};

test_fun = @(x) sin(x);                 element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[cos(x_val) 0 0 0 0]},5};
test_fun = @(x) sinh(x);                element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[cosh(x_val) 0 0 0 0]},5};
test_fun = @(x) asin(x);                element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[1/sqrt(1-x_val^2) 0 0 0 0]},5};
test_fun = @(x) asinh(x);               element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[1/sqrt(x_val^2+1) 0 0 0 0]},5};

test_fun = @(x) cos(x);                 element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[-sin(x_val) 0 0 0 0]},5};
test_fun = @(x) cosh(x);                element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[sinh(x_val) 0 0 0 0]},5};
test_fun = @(x) acos(x);                element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[-1/sqrt(1-x_val^2) 0 0 0 0]},5};
test_fun = @(x) acosh(x);               element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[ 1/sqrt(x_val^2-1) 0 0 0 0]},5};

test_fun = @(x) tan(x);                 element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[sec(x_val)^2 0 0 0 0]},5};
test_fun = @(x) tanh(x);                element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[1-tanh(x_val)^2 0 0 0 0]},5};
test_fun = @(x) atan(x);                element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[1/(1+x_val^2) 0 0 0 0]},5};
test_fun = @(x) atanh(x);               element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[1/(1-x_val^2) 0 0 0 0]},5};

test_fun = @(x) cot(x);                 element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[-csc(x_val)^2 0 0 0 0]},5};
test_fun = @(x) coth(x);                element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[1 - coth(x_val)^2 0 0 0 0]},5};
test_fun = @(x) acot(x);                element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[-1/(1+x_val^2) 0 0 0 0]},5};
test_fun = @(x) acoth(x);               element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[1/(1-x_val^2) 0 0 0 0]},5};

test_fun = @(x) csc(x);                 element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[-cot(x_val)*csc(x_val) 0 0 0 0]},5};
test_fun = @(x) csch(x);                element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[-coth(x_val)*csch(x_val) 0 0 0 0]},5};
test_fun = @(x) acsc(x);                element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[-1/(abs(x_val)*sqrt(x_val^2-1)) 0 0 0 0]},5};
test_fun = @(x) acsch(x);               element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[-1/(abs(x_val)*sqrt(1+x_val^2)) 0 0 0 0]},5};

test_fun = @(x) sec(x);                 element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[sec(x_val)*tan(x_val) 0 0 0 0]},5};
test_fun = @(x) sech(x);                element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[-tanh(x_val)*sech(x_val) 0 0 0 0]},5};
test_fun = @(x) asec(x);                element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[1/(abs(x_val)*sqrt(x_val^2-1)) 0 0 0 0]},5};
test_fun = @(x) asech(x);               element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[-1/(abs(x_val)*sqrt(1-x_val^2)) 0 0 0 0]},5};

test_fun = @(x) log(x);                 element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[1/x_val 0 0 0 0]},5};
test_fun = @(x) gamma(x);               element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[gamma(x_val)*psi(x_val) 0 0 0 0]},5};
test_fun = @(x) x^3.41;                 element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[3.41*x_val^2.41 0 0 0 0]},5};
test_fun = @(x) x^0;                    element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[0 0 0 0 0]},5};
test_fun = @(x) x^(-1);                 element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[-1/(x_val^2) 0 0 0 0]},5};

test_fun = @(x) 3.456^x;                element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[log(3.456)*3.456^x_val 0 0 0 0]},5};
test_fun = @(x) x.^3.41;                element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[3.41*x_val^2.41 0 0 0 0]},5};
test_fun = @(x) x.^0;                   element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[0 0 0 0 0]},5};
test_fun = @(x) x.^(-1);                element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[-1/(x_val^2) 0 0 0 0]},5};
test_fun = @(x) 3.456.^x;               element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[log(3.456)*3.456^x_val 0 0 0 0]},5};

test_fun = @(x) sqrt(x);                element1{length(element1)+1}       = {test_fun,{x},test_fun(x_val)          ,{[0.5/sqrt(x_val) 0 0 0 0]},5};

%% MATRIX ARITHMETIC OPERATION - SECOND TEST SERIES
%function to test                   Save in as one element                 func      in     corr val                     corr derivative
test_fun = @(x, y, z, w, q) 2.4*x+4.3*y+1.5*z+5.6*w+9.1*q;                 element2{1}                          = {test_fun,{x y z w q}, test_fun(x_val,y_val,z_val,w_val,q_val), {[2.4 4.3 1.5 5.6 9.1]},5};
test_fun = @(x, y, z, w, q) -5.1*x-1.2*y-3.5*z-7.4*w-7.3*q;                element2{length(element2)+1}         = {test_fun,{x y z w q}, test_fun(x_val,y_val,z_val,w_val,q_val) ,{[-5.1 -1.2 -3.5 -7.4 -7.3]},5};
test_fun = @(x, y, z, w, q) 2.4*x-1.2*y-3.5*z+5.6*w-7.3*q;                 element2{length(element2)+1}         = {test_fun,{x y z w q}, test_fun(x_val,y_val,z_val,w_val,q_val) ,{[2.4 -1.2 -3.5 5.6 -7.3]},5};
test_fun = @(x, y, z, w, q) (2.4*x-1.2*y)-3.5*z+(5.6*w-7.3*q);             element2{length(element2)+1}         = {test_fun,{x y z w q}, test_fun(x_val,y_val,z_val,w_val,q_val) ,{[2.4 -1.2 -3.5 5.6 -7.3]},5};
test_fun = @(x, y, z, w, q) (-5.1*x-1.2*y)-(3.5*z-7.4*w-7.3*q);            element2{length(element2)+1}         = {test_fun,{x y z w q}, test_fun(x_val,y_val,z_val,w_val,q_val) ,{[-5.1 -1.2 -3.5 +7.4 +7.3]},5};

test_fun = @(x,y) x*y;                  element2{length(element2)+1}       = {test_fun,{x y}, test_fun(x_val,y_val) ,{[y_val x_val 0 0 0]},5};
test_fun = @(x,y) x.*y;                 element2{length(element2)+1}       = {test_fun,{x y}, test_fun(x_val,y_val) ,{[y_val x_val 0 0 0]},5};
test_fun = @(x,y) x/y;                  element2{length(element2)+1}       = {test_fun,{x y}, test_fun(x_val,y_val) ,{[1/y_val -x_val/(y_val^2) 0 0 0]},5};
test_fun = @(x,y) x./y;                 element2{length(element2)+1}       = {test_fun,{x y}, test_fun(x_val,y_val) ,{[1/y_val -x_val/(y_val^2) 0 0 0]},5};
test_fun = @(x,y) x./y;                 element2{length(element2)+1}       = {test_fun,{x y}, test_fun(x_val,y_val) ,{[1/y_val -x_val/(y_val^2) 0 0 0]},5};

test_fun = @(x,y) x^y;                  element2{length(element2)+1}       = {test_fun,{x y}, test_fun(x_val,y_val) ,{[y_val*x_val^(y_val-1) log(x_val)*x_val^y_val 0 0 0]},5};
test_fun = @(x,y) x.^y;                 element2{length(element2)+1}       = {test_fun,{x y}, test_fun(x_val,y_val) ,{[y_val*x_val^(y_val-1) log(x_val)*x_val^y_val 0 0 0]},5};
test_fun = @(x,y,z) x + y*z;            element2{length(element2)+1}       = {test_fun,{x y z}, test_fun(x_val, y_val, z_val)         ,{[1 z_val y_val 0 0]},5};
test_fun = @(x) sin(cos(exp(tan(x))));  element2{length(element2)+1}       = {test_fun,{x}, test_fun(x_val)         ,{[-exp(tan(x_val)) * sec(x_val)^2 * sin(exp(tan(x_val))) *cos(cos(exp(tan(x_val)))) 0 0 0 0 ]},5};
test_fun = @(x) 1.14*x;                 element2{length(element2)+1}       = {test_fun,{A}, test_fun(A_val)         ,{[1.14 0 0 0] [0 0 1.14 0];[0 1.14 0 0] [0 0 0 1.14]},4};
test_fun = @(x,y) x + y;                element2{length(element2)+1}       = {test_fun,{x j}, test_fun(x_val,j_val)         ,{[1 0 0 0 0]},5};

%% MATRIX ARITHMETIC OPERATION - THIRD TEST SERIES

test_fun = @(x) x+2.41;                 element3{1}                        = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix((ones(size(B_val)))),n_in_B};
test_fun = @(x) x-2.81;                 element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix((ones(size(B_val)))),n_in_B};
test_fun = @(x) 6.41+x;                 element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix((ones(size(B_val)))),n_in_B};
test_fun = @(x) 6.5 -x;                 element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix((-ones(size(B_val)))),n_in_B};
test_fun = @(x) -x;                     element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix((-ones(size(B_val)))),n_in_B};

test_fun = @(x) x*9.2;                  element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix((9.2*ones(size(B_val)))),n_in_B};
test_fun = @(x) x.*2.11;                element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix((2.11*ones(size(B_val)))),n_in_B};
test_fun = @(x) 4.2*x;                  element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix((4.2*ones(size(B_val)))),n_in_B};
test_fun = @(x) 2.46.*x;                element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix((2.46*ones(size(B_val)))),n_in_B};

test_fun = @(x) x/0.4;                  element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix((1/0.4*ones(size(B_val)))),n_in_B};
test_fun = @(x) x./2.8;                 element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix((1/2.8*ones(size(B_val)))),n_in_B};
test_fun = @(x) 5.46./x;                element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(-5.46./(B_val.^2)),n_in_B};

test_fun = @(x) x.^4.51;                element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(4.51.*B_val.^3.51),n_in_B};
test_fun = @(x) x.^0;                   element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(zeros(size(B_val))),n_in_B};
test_fun = @(x) x.^-2.71;               element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(-2.71.*B_val.^(-3.71)),n_in_B};

test_fun = @(x) uminus(x);              element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(-(ones(size(B_val)))),n_in_B};
test_fun = @(x) uminus(x);              element3{length(element3)+1}       = {test_fun,{b}, test_fun(b_val)         ,der_from_matrix(-(ones(size(b_val)))),n_in_b};
test_fun = @(x) uminus(x);              element3{length(element3)+1}       = {test_fun,{b'},test_fun(b_val')        ,der_from_matrix(-(ones(size(b_val')))),n_in_b};

test_fun = @(x) exp(x);                 element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(exp(B_val)),n_in_B};
test_fun = @(x) exp(x);                 element3{length(element3)+1}       = {test_fun,{b}, test_fun(b_val)         ,der_from_matrix(exp(b_val)),n_in_b};
test_fun = @(x) exp(x);                 element3{length(element3)+1}       = {test_fun,{b'},test_fun(b_val')        ,der_from_matrix(exp(b_val')),n_in_b};

test_fun = @(x) log(x);                 element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(1./B_val),n_in_B};
test_fun = @(x) log(x);                 element3{length(element3)+1}       = {test_fun,{b}, test_fun(b_val)         ,der_from_matrix(1./b_val),n_in_b};
test_fun = @(x) log(x);                 element3{length(element3)+1}       = {test_fun,{b'},test_fun(b_val')        ,der_from_matrix(1./b_val'),n_in_b};

test_fun = @(x) sqrt(x);                element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(1./(2.*sqrt(B_val))),n_in_B};
test_fun = @(x) sqrt(x);                element3{length(element3)+1}       = {test_fun,{b}, test_fun(b_val)         ,der_from_matrix(1./(2.*sqrt(b_val))),n_in_b};
test_fun = @(x) sqrt(x);                element3{length(element3)+1}       = {test_fun,{b'},test_fun(b_val')        ,der_from_matrix(1./(2.*sqrt(b_val'))),n_in_b};

test_fun = @(x) sin(x);                 element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(cos(B_val)),n_in_B};
test_fun = @(x) sin(x);                 element3{length(element3)+1}       = {test_fun,{b}, test_fun(b_val)         ,der_from_matrix(cos(b_val)),n_in_b};
test_fun = @(x) sin(x);                 element3{length(element3)+1}       = {test_fun,{b'},test_fun(b_val')        ,der_from_matrix(cos(b_val')),n_in_b};
test_fun = @(x) sinh(x);                element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(cosh(B_val)),n_in_B};
test_fun = @(x) sinh(x);                element3{length(element3)+1}       = {test_fun,{b}, test_fun(b_val)         ,der_from_matrix(cosh(b_val)),n_in_b};
test_fun = @(x) sinh(x);                element3{length(element3)+1}       = {test_fun,{b'},test_fun(b_val')        ,der_from_matrix(cosh(b_val')),n_in_b};
test_fun = @(x) asin(x);                element3{length(element3)+1}       = {test_fun,{Bs}, test_fun(Bs_val)       ,der_from_matrix(1./sqrt(1-Bs_val.^2)) ,n_in_B}; 
test_fun = @(x) asin(x);                element3{length(element3)+1}       = {test_fun,{bs}, test_fun(bs_val)       ,der_from_matrix(1./sqrt(1-bs_val.^2)) ,n_in_b};
test_fun = @(x) asin(x);                element3{length(element3)+1}       = {test_fun,{bs'},test_fun(bs_val')      ,der_from_matrix(1./sqrt(1-bs_val'.^2)),n_in_b};
test_fun = @(x) asinh(x);               element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(1./sqrt(B_val.^2+1)),n_in_B};
test_fun = @(x) asinh(x);               element3{length(element3)+1}       = {test_fun,{b}, test_fun(b_val)         ,der_from_matrix(1./sqrt(b_val.^2+1)),n_in_b};
test_fun = @(x) asinh(x);               element3{length(element3)+1}       = {test_fun,{b'},test_fun(b_val')        ,der_from_matrix(1./sqrt(b_val'.^2+1)),n_in_b};

test_fun = @(x) cos(x);                 element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(-sin(B_val)),n_in_B};
test_fun = @(x) cos(x);                 element3{length(element3)+1}       = {test_fun,{b}, test_fun(b_val)         ,der_from_matrix(-sin(b_val)),n_in_b};
test_fun = @(x) cos(x);                 element3{length(element3)+1}       = {test_fun,{b'},test_fun(b_val')        ,der_from_matrix(-sin(b_val')),n_in_b};
test_fun = @(x) cosh(x);                element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(sinh(B_val)),n_in_B};
test_fun = @(x) cosh(x);                element3{length(element3)+1}       = {test_fun,{b}, test_fun(b_val)         ,der_from_matrix(sinh(b_val)),n_in_b};
test_fun = @(x) cosh(x);                element3{length(element3)+1}       = {test_fun,{b'},test_fun(b_val')        ,der_from_matrix(sinh(b_val')),n_in_b};
test_fun = @(x) acos(x);                element3{length(element3)+1}       = {test_fun,{Bs}, test_fun(Bs_val)       ,der_from_matrix(-1./sqrt(1-Bs_val.^2)),n_in_B};
test_fun = @(x) acos(x);                element3{length(element3)+1}       = {test_fun,{bs}, test_fun(bs_val)       ,der_from_matrix(-1./sqrt(1-bs_val.^2)),n_in_b};
test_fun = @(x) acos(x);                element3{length(element3)+1}       = {test_fun,{bs'},test_fun(bs_val')      ,der_from_matrix(-1./sqrt(1-bs_val'.^2)),n_in_b};
test_fun = @(x) acosh(x);               element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix( 1./sqrt(B_val.^2-1) ),n_in_B};
test_fun = @(x) acosh(x);               element3{length(element3)+1}       = {test_fun,{b}, test_fun(b_val)         ,der_from_matrix( 1./sqrt(b_val.^2-1) ),n_in_b};
test_fun = @(x) acosh(x);               element3{length(element3)+1}       = {test_fun,{b'},test_fun(b_val')        ,der_from_matrix( 1./sqrt(b_val'.^2-1) ),n_in_b};

test_fun = @(x) tan(x);                 element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(sec(B_val).^2)    ,n_in_B};
test_fun = @(x) tan(x);                 element3{length(element3)+1}       = {test_fun,{b}, test_fun(b_val)         ,der_from_matrix(sec(b_val).^2)    ,n_in_b};
test_fun = @(x) tan(x);                 element3{length(element3)+1}       = {test_fun,{b'},test_fun(b_val')        ,der_from_matrix(sec(b_val').^2)   ,n_in_b};
test_fun = @(x) tanh(x);                element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(1-tanh(B_val).^2) ,n_in_B};
test_fun = @(x) tanh(x);                element3{length(element3)+1}       = {test_fun,{b}, test_fun(b_val)         ,der_from_matrix(1-tanh(b_val).^2) ,n_in_b};
test_fun = @(x) tanh(x);                element3{length(element3)+1}       = {test_fun,{b'},test_fun(b_val')        ,der_from_matrix(1-tanh(b_val').^2),n_in_b};
test_fun = @(x) atan(x);                element3{length(element3)+1}       = {test_fun,{Bs}, test_fun(Bs_val)       ,der_from_matrix(1./(1+Bs_val.^2))   ,n_in_B};
test_fun = @(x) atan(x);                element3{length(element3)+1}       = {test_fun,{bs}, test_fun(bs_val)       ,der_from_matrix(1./(1+bs_val.^2))   ,n_in_b};
test_fun = @(x) atan(x);                element3{length(element3)+1}       = {test_fun,{bs'},test_fun(bs_val')      ,der_from_matrix(1./(1+bs_val'.^2))  ,n_in_b};
test_fun = @(x) atanh(x);               element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(1./(1-B_val.^2) ) ,n_in_B};
test_fun = @(x) atanh(x);               element3{length(element3)+1}       = {test_fun,{b}, test_fun(b_val)         ,der_from_matrix(1./(1-b_val.^2) ) ,n_in_b};
test_fun = @(x) atanh(x);               element3{length(element3)+1}       = {test_fun,{b'},test_fun(b_val')        ,der_from_matrix(1./(1-b_val'.^2) ),n_in_b};

test_fun = @(x) sec(x);                 element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(sec(B_val).*tan(B_val))  ,n_in_B};
test_fun = @(x) sec(x);                 element3{length(element3)+1}       = {test_fun,{b}, test_fun(b_val)         ,der_from_matrix(sec(b_val).*tan(b_val))  ,n_in_b};
test_fun = @(x) sec(x);                 element3{length(element3)+1}       = {test_fun,{b'},test_fun(b_val')        ,der_from_matrix(sec(b_val').*tan(b_val')),n_in_b};
test_fun = @(x) sech(x);                element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(-tanh(B_val).*sech(B_val)),n_in_B};
test_fun = @(x) sech(x);                element3{length(element3)+1}       = {test_fun,{b}, test_fun(b_val)         ,der_from_matrix(-tanh(b_val).*sech(b_val)),n_in_b};
test_fun = @(x) sech(x);                element3{length(element3)+1}       = {test_fun,{b'},test_fun(b_val')        ,der_from_matrix(-tanh(b_val').*sech(b_val')),n_in_b};
test_fun = @(x) asec(x);                element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(1./(abs(B_val).*sqrt(B_val.^2-1))),n_in_B};
test_fun = @(x) asec(x);                element3{length(element3)+1}       = {test_fun,{b}, test_fun(b_val)         ,der_from_matrix(1./(abs(b_val).*sqrt(b_val.^2-1))),n_in_b};
test_fun = @(x) asec(x);                element3{length(element3)+1}       = {test_fun,{b'},test_fun(b_val')        ,der_from_matrix(1./(abs(b_val').*sqrt(b_val'.^2-1))),n_in_b};
test_fun = @(x) asech(x);               element3{length(element3)+1}       = {test_fun,{Bs}, test_fun(Bs_val)       ,der_from_matrix(-1./(Bs_val.*sqrt(1-Bs_val.^2)) ),n_in_B};
test_fun = @(x) asech(x);               element3{length(element3)+1}       = {test_fun,{bs}, test_fun(bs_val)       ,der_from_matrix(-1./(bs_val.*sqrt(1-bs_val.^2)) ),n_in_b};
test_fun = @(x) asech(x);               element3{length(element3)+1}       = {test_fun,{bs'},test_fun(bs_val')      ,der_from_matrix(-1./(bs_val'.*sqrt(1-bs_val'.^2)) ),n_in_b};

test_fun = @(x) cot(x);                 element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(-csc(B_val).^2)   ,n_in_B};
test_fun = @(x) cot(x);                 element3{length(element3)+1}       = {test_fun,{b}, test_fun(b_val)         ,der_from_matrix(-csc(b_val).^2)   ,n_in_b};
test_fun = @(x) cot(x);                 element3{length(element3)+1}       = {test_fun,{b'},test_fun(b_val')        ,der_from_matrix(-csc(b_val').^2)  ,n_in_b};
test_fun = @(x) coth(x);                element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(1-coth(B_val).^2) ,n_in_B};
test_fun = @(x) coth(x);                element3{length(element3)+1}       = {test_fun,{b}, test_fun(b_val)         ,der_from_matrix(1-coth(b_val).^2) ,n_in_b};
test_fun = @(x) coth(x);                element3{length(element3)+1}       = {test_fun,{b'},test_fun(b_val')        ,der_from_matrix(1-coth(b_val').^2),n_in_b};
test_fun = @(x) acot(x);                element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(-1./(1+B_val.^2))   ,n_in_B};
test_fun = @(x) acot(x);                element3{length(element3)+1}       = {test_fun,{b}, test_fun(b_val)         ,der_from_matrix(-1./(1+b_val.^2))   ,n_in_b};
test_fun = @(x) acot(x);                element3{length(element3)+1}       = {test_fun,{b'},test_fun(b_val')        ,der_from_matrix(-1./(1+b_val'.^2))  ,n_in_b};
test_fun = @(x) acoth(x);               element3{length(element3)+1}       = {test_fun,{Bs}, test_fun(Bs_val)       ,der_from_matrix(1./(1-Bs_val.^2)) ,n_in_B};
test_fun = @(x) acoth(x);               element3{length(element3)+1}       = {test_fun,{bs}, test_fun(bs_val)       ,der_from_matrix(1./(1-bs_val.^2)) ,n_in_b};
test_fun = @(x) acoth(x);               element3{length(element3)+1}       = {test_fun,{bs'},test_fun(bs_val')      ,der_from_matrix(1./(1-bs_val'.^2)),n_in_b};

test_fun = @(x) csc(x);                 element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(-cot(B_val).*csc(B_val) )  ,n_in_B};
test_fun = @(x) csc(x);                 element3{length(element3)+1}       = {test_fun,{b}, test_fun(b_val)         ,der_from_matrix(-cot(b_val).*csc(b_val) )  ,n_in_b};
test_fun = @(x) csc(x);                 element3{length(element3)+1}       = {test_fun,{b'},test_fun(b_val')        ,der_from_matrix(-cot(b_val').*csc(b_val')) ,n_in_b};
test_fun = @(x) csch(x);                element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(-coth(B_val).*csch(B_val) )  ,n_in_B};
test_fun = @(x) csch(x);                element3{length(element3)+1}       = {test_fun,{b}, test_fun(b_val)         ,der_from_matrix(-coth(b_val).*csch(b_val) )  ,n_in_b};
test_fun = @(x) csch(x);                element3{length(element3)+1}       = {test_fun,{b'},test_fun(b_val')        ,der_from_matrix(-coth(b_val').*csch(b_val')) ,n_in_b};
test_fun = @(x) acsc(x);                element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(-1./(abs(B_val).*sqrt(B_val.^2-1))  )  ,n_in_B};
test_fun = @(x) acsc(x);                element3{length(element3)+1}       = {test_fun,{b}, test_fun(b_val)         ,der_from_matrix(-1./(abs(b_val).*sqrt(b_val.^2-1))  )  ,n_in_b};
test_fun = @(x) acsc(x);                element3{length(element3)+1}       = {test_fun,{b'},test_fun(b_val')        ,der_from_matrix(-1./(abs(b_val').*sqrt(b_val'.^2-1))  ) ,n_in_b};
test_fun = @(x) acsch(x);               element3{length(element3)+1}       = {test_fun,{B}, test_fun(B_val)         ,der_from_matrix(-1./(abs(B_val).*sqrt(1+B_val.^2)) )  ,n_in_B};
test_fun = @(x) acsch(x);               element3{length(element3)+1}       = {test_fun,{b}, test_fun(b_val)         ,der_from_matrix(-1./(abs(b_val).*sqrt(1+b_val.^2)) )  ,n_in_b};
test_fun = @(x) acsch(x);               element3{length(element3)+1}       = {test_fun,{b'},test_fun(b_val')        ,der_from_matrix(-1./(abs(b_val').*sqrt(1+b_val'.^2)) ) ,n_in_b};

%% MATRIX ARITHMETIC OPERATION - FOURTH TEST SERIES
test_fun = @(x) [2 7 6;8 6 3].*x;       element4{1}                        = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix([2 7 6;8 6 3]),n_in_B};
test_fun = @(x) x.*[2 7 6;8 6 3];       element4{length(element4)+1}       = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix([2 7 6;8 6 3]),n_in_B};

test_fun = @(x) x./[2 7 6;8 6 3];       element4{length(element4)+1}       = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix(1./[2 7 6;8 6 3]),n_in_B};
test_fun = @(x) [2 7 6;8 6 3]./x;       element4{length(element4)+1}       = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix(-[2 7 6;8 6 3]./(B_val.^2)),n_in_B};

test_fun = @(x) x + [2 7 6;8 6 3];      element4{length(element4)+1}       = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix(ones(size(B_val))),n_in_B};
test_fun = @(x) [2 7 6;8 6 3] + x;      element4{length(element4)+1}       = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix(ones(size(B_val))),n_in_B};

test_fun = @(x) x - [2 7 6;8 6 3];      element4{length(element4)+1}       = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix(ones(size(B_val))),n_in_B};
test_fun = @(x) [2 7 6;8 6 3] - x;      element4{length(element4)+1}       = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix(-ones(size(B_val))),n_in_B};

test_fun = @(x) x.^[2 7 6;8 6 3];       element4{length(element4)+1}       = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix([2 7 6;8 6 3].*B_val.^([2 7 6;8 6 3]-1)),n_in_B};
test_fun = @(x) [2 7 6;8 6 3].^x;       element4{length(element4)+1}       = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix(log([2 7 6;8 6 3]).*[2 7 6;8 6 3].^(B_val)),n_in_B};

test_fun = @(x) x*[4; 3; 7];            element4{length(element4)+1}       = {test_fun,{B}, test_fun(B_val)          ,{[4 0 3 0 7 0];[0 4 0 3 0 7]},n_in_B};
test_fun = @(x) [3 5]*x;                element4{length(element4)+1}       = {test_fun,{B}, test_fun(B_val)          ,{[3 5 0 0 0 0],[0 0 3 5 0 0],[0 0 0 0 3 5]},n_in_B};

test_fun = @(x) [3 5].*x;               element4{length(element4)+1}       = {test_fun,{x}, test_fun(x_val)          ,{[3 0 0 0 0],[5 0 0 0 0]},5};

%% MATRIX ARITHMETIC OPERATION - FIFTH TEST SERIES
test_fun = @(x) x.*x;                   element5{1}                        = {test_fun,{A}, test_fun([A_val])        ,der_from_matrix(A_val.*2),n_in_A};
test_fun = @(x) x*x;                    element5{length(element5)+1}       = {test_fun,{A}, test_fun([A_val])        ,{[6 2 4 0],[2 0 8 2];[4 8 0 4],[0 2 4 10]},n_in_A};
test_fun = @(x,y) x + y;                element5{length(element5)+1}       = {test_fun,{B C},test_fun(B_val, C_val)  ,{[1 0 0 0 0 0 1 0 0 0 0 0],[0 0 1 0 0 0 0 0 1 0 0 0],[0 0 0 0 1 0 0 0 0 0 1 0];[0 1 0 0 0 0 0 1 0 0 0 0],[0 0 0 1 0 0 0 0 0 1 0 0],[0 0 0 0 0 1 0 0 0 0 0 1]},n_in_B*2};
test_fun = @(x,y) x - y;                element5{length(element5)+1}       = {test_fun,{B C},test_fun(B_val, C_val)  ,{[1 0 0 0 0 0 -1 0 0 0 0 0],[0 0 1 0 0 0 0 0 -1 0 0 0],[0 0 0 0 1 0 0 0 0 0 -1 0];[0 1 0 0 0 0 0 -1 0 0 0 0],[0 0 0 1 0 0 0 0 0 -1 0 0],[0 0 0 0 0 1 0 0 0 0 0 -1]},n_in_B*2};
test_fun = @(x,y) x./ y;                element5{length(element5)+1}       = {test_fun,{C B},test_fun(C_val, B_val)  ,{[-C_val(1,1)/(B_val(1,1).^2) 0 0 0 0 0 1/B_val(1,1) 0 0 0 0 0],...
                                                                                                                       [0 0 -C_val(1,2)/(B_val(1,2).^2) 0 0 0 0 0 1/B_val(1,2) 0 0 0],...
                                                                                                                       [0 0 0 0 -C_val(1,3)/(B_val(1,3).^2) 0 0 0 0 0 1/B_val(1,3) 0];...
                                                                                                                       [0 -C_val(2,1)/(B_val(2,1).^2) 0 0 0 0 0 1/B_val(2,1) 0 0 0 0],...
                                                                                                                       [0 0 0 -C_val(2,2)/(B_val(2,2).^2) 0 0 0 0 0 1/B_val(2,2) 0 0],...
                                                                                                                       [0 0 0 0 0 -C_val(2,3)/(B_val(2,3).^2) 0 0 0 0 0 1/B_val(2,3)]},n_in_B*2};

test_fun = @(x,y) x.^y;                 element5{length(element5)+1}       = {test_fun,{B C},[test_fun(B_val, C_val)],{[C_val(1,1).*(B_val(1,1).^(C_val(1,1)-1)) 0 0 0 0 0 log(B_val(1,1)).*B_val(1,1).^(C_val(1,1)) 0 0 0 0 0],...
                                                                                                                       [0 0 C_val(1,2).*(B_val(1,2).^(C_val(1,2)-1)) 0 0 0 0 0 log(B_val(1,2)).*B_val(1,2).^(C_val(1,2)) 0 0 0],...
                                                                                                                       [0 0 0 0 C_val(1,3).*(B_val(1,3).^(C_val(1,3)-1)) 0 0 0 0 0 log(B_val(1,3)).*B_val(1,3).^(C_val(1,3)) 0];...
                                                                                                                       [0 C_val(2,1).*(B_val(2,1).^(C_val(2,1)-1)) 0 0 0 0 0 log(B_val(2,1)).*B_val(2,1).^(C_val(2,1)) 0 0 0 0],...
                                                                                                                       [0 0 0 C_val(2,2).*(B_val(2,2).^(C_val(2,2)-1)) 0 0 0 0 0 log(B_val(2,2)).*B_val(2,2).^(C_val(2,2)) 0 0],...
                                                                                                                       [0 0 0 0 0 C_val(2,3).*(B_val(2,3).^(C_val(2,3)-1)) 0 0 0 0 0 log(B_val(2,3)).*B_val(2,3).^(C_val(2,3))]},n_in_B*2};                                                                                                                      
                                                                                                                  
test_fun = @(x) x*Fb;                   element5{length(element5)+1}       = {test_fun,[F], [test_fun([F_val])]       ,M_v_der_out(F_val,Fb),n_in_F};
test_fun = @(x) Fb'*x;                  element5{length(element5)+1}       = {test_fun,[F], [test_fun([F_val])]       ,v_M_der_out(F_val,Fb'),n_in_F};
test_fun = @(x) FAb*x;                  element5{length(element5)+1}       = {test_fun,[F], [test_fun([F_val])]       ,M_M_der_out(FAb,F_val,'first'),n_in_F};
test_fun = @(x) x*FAa;                  element5{length(element5)+1}       = {test_fun,[F], [test_fun([F_val])]       ,M_M_der_out(F_val,FAa,'second'),n_in_F};
test_fun = @(x,y) x*y;                  element5{length(element5)+1}       = {test_fun,{G H}, [test_fun(G_val, H_val)],M_M_der_out(G_val,H_val,'none'),n_in_G+n_in_H};
test_fun = @(x) x + [8 9 10];           element5{length(element5)+1}       = {test_fun,{x}, [test_fun(x_val)],{[1 0 0 0 0],[1 0 0 0 0],[1 0 0 0 0]},5};
test_fun = @(x) [8 9 10] + x;           element5{length(element5)+1}       = {test_fun,{x}, [test_fun(x_val)],{[1 0 0 0 0],[1 0 0 0 0],[1 0 0 0 0]},5};
test_fun = @(x,y) x.*y;                 element5{length(element5)+1}       = {test_fun,{x,yy}, [test_fun(x_val,yy_val)],{[yy_val(1) x_val 0 0 0],[yy_val(2) 0 x_val 0 0]},5};

%% MATRIX ARITHMETIC OPERATION - SIXTH TEST SERIES PLUS(+)
test_fun = @(x) x + [2 7 6];            element6{1}                        = {test_fun,{x}, test_fun(x_val)          ,{[1 0 0 0 0],[1 0 0 0 0],[1 0 0 0 0]},5};
test_fun = @(x) x + [2;7;6];            element6{length(element6)+1}       = {test_fun,{x}, test_fun(x_val)          ,{[1 0 0 0 0];[1 0 0 0 0];[1 0 0 0 0]},5};
test_fun = @(x) x + [2 7 6;8 6 3];      element6{length(element6)+1}       = {test_fun,{x}, test_fun(x_val)          ,{[1 0 0 0 0],[1 0 0 0 0],[1 0 0 0 0];[1 0 0 0 0],[1 0 0 0 0],[1 0 0 0 0]},5};

test_fun = @(x) x + 8;                  element6{length(element6)+1}       = {test_fun,{b}, test_fun(b_val)          ,{[1 0 0 0 0 0],[0 1 0 0 0 0],[0 0 1 0 0 0],[0 0 0 1 0 0],[0 0 0 0 1 0],[0 0 0 0 0 1]},6};
test_fun = @(x) x + [2 7 6 8 6 3];      element6{length(element6)+1}       = {test_fun,{b}, test_fun(b_val)          ,{[1 0 0 0 0 0],[0 1 0 0 0 0],[0 0 1 0 0 0],[0 0 0 1 0 0],[0 0 0 0 1 0],[0 0 0 0 0 1]},6};

test_fun = @(x) x + 8;                  element6{length(element6)+1}       = {test_fun,{b'}, test_fun(b_val')        ,{[1 0 0 0 0 0];[0 1 0 0 0 0];[0 0 1 0 0 0];[0 0 0 1 0 0];[0 0 0 0 1 0];[0 0 0 0 0 1]},6};
test_fun = @(x) x + [2 7 6 8 6 3]';     element6{length(element6)+1}       = {test_fun,{b'}, test_fun(b_val')        ,{[1 0 0 0 0 0];[0 1 0 0 0 0];[0 0 1 0 0 0];[0 0 0 1 0 0];[0 0 0 0 1 0];[0 0 0 0 0 1]},6};

test_fun = @(x) x+6.41;                 element6{length(element6)+1}       = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix((ones(size(B_val)))),n_in_B};
test_fun = @(x) x + [2 7 6;8 6 3];      element6{length(element6)+1}       = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix(ones(size(B_val))),n_in_B};

%

test_fun = @(x) [2 7 6] + x;            element6{length(element6)+1}       = {test_fun,{x}, test_fun(x_val)          ,{[1 0 0 0 0],[1 0 0 0 0],[1 0 0 0 0]},5};
test_fun = @(x) [2;7;6] + x;            element6{length(element6)+1}       = {test_fun,{x}, test_fun(x_val)          ,{[1 0 0 0 0];[1 0 0 0 0];[1 0 0 0 0]},5};
test_fun = @(x) [2 7 6;8 6 3] + x;      element6{length(element6)+1}       = {test_fun,{x}, test_fun(x_val)          ,{[1 0 0 0 0],[1 0 0 0 0],[1 0 0 0 0];[1 0 0 0 0],[1 0 0 0 0],[1 0 0 0 0]},5};

test_fun = @(x) 8 + x;                  element6{length(element6)+1}       = {test_fun,{b}, test_fun(b_val)          ,{[1 0 0 0 0 0],[0 1 0 0 0 0],[0 0 1 0 0 0],[0 0 0 1 0 0],[0 0 0 0 1 0],[0 0 0 0 0 1]},6};
test_fun = @(x) [2 7 6 8 6 3] + x;      element6{length(element6)+1}       = {test_fun,{b}, test_fun(b_val)          ,{[1 0 0 0 0 0],[0 1 0 0 0 0],[0 0 1 0 0 0],[0 0 0 1 0 0],[0 0 0 0 1 0],[0 0 0 0 0 1]},6};

test_fun = @(x) 8 + x;                  element6{length(element6)+1}       = {test_fun,{b'}, test_fun(b_val')        ,{[1 0 0 0 0 0];[0 1 0 0 0 0];[0 0 1 0 0 0];[0 0 0 1 0 0];[0 0 0 0 1 0];[0 0 0 0 0 1]},6};
test_fun = @(x) [2 7 6 8 6 3]' + x;     element6{length(element6)+1}       = {test_fun,{b'}, test_fun(b_val')        ,{[1 0 0 0 0 0];[0 1 0 0 0 0];[0 0 1 0 0 0];[0 0 0 1 0 0];[0 0 0 0 1 0];[0 0 0 0 0 1]},6};

test_fun = @(x) 8 + x;                  element6{length(element6)+1}       = {test_fun,{b'}, test_fun(b_val')        ,{[1 0 0 0 0 0];[0 1 0 0 0 0];[0 0 1 0 0 0];[0 0 0 1 0 0];[0 0 0 0 1 0];[0 0 0 0 0 1]},6};
test_fun = @(x) [2 7 6 8 6 3]' + x;     element6{length(element6)+1}       = {test_fun,{b'}, test_fun(b_val')        ,{[1 0 0 0 0 0];[0 1 0 0 0 0];[0 0 1 0 0 0];[0 0 0 1 0 0];[0 0 0 0 1 0];[0 0 0 0 0 1]},6};

test_fun = @(x) 6.41 + x;               element6{length(element6)+1}       = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix((ones(size(B_val)))),n_in_B};
test_fun = @(x) [2 7 6;8 6 3] + x;      element6{length(element6)+1}       = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix(ones(size(B_val))),n_in_B};

%

test_fun = @(x,y) x + y;                element6{length(element6)+1}       = {test_fun,{x,y},   test_fun(x_val,y_val)   ,{[1 1 0 0 0]},5};
test_fun = @(x,y) x + y;                element6{length(element6)+1}       = {test_fun,{x,yy},  test_fun(x_val,yy_val)  ,{[1 1 0 0 0],[1 0 1 0 0]},5};
test_fun = @(x,y) x + y;                element6{length(element6)+1}       = {test_fun,{x,yy'}, test_fun(x_val,yy_val') ,{[1 1 0 0 0];[1 0 1 0 0]},5};
test_fun = @(x,y) x + y;                element6{length(element6)+1}       = {test_fun,{x,B},   test_fun(x_val,B_val)   ,{[2 0 0 0 0 0],[1 0 1 0 0 0],[1 0 0 0 1 0];[1 1 0 0 0 0],[1 0 0 1 0 0],[1 0 0 0 0 1]},6};

test_fun = @(x,y) x + y;                element6{length(element6)+1}       = {test_fun,{yy,x},  test_fun(yy_val,x_val)  ,{[1 1 0 0 0],[1 0 1 0 0]},5};
test_fun = @(x,y) x + y;                element6{length(element6)+1}       = {test_fun,{yy',x}, test_fun(yy_val',x_val) ,{[1 1 0 0 0];[1 0 1 0 0]},5};
test_fun = @(x,y) x + y;                element6{length(element6)+1}       = {test_fun,{B,x},   test_fun(B_val,x_val)   ,{[2 0 0 0 0 0],[1 0 1 0 0 0],[1 0 0 0 1 0];[1 1 0 0 0 0],[1 0 0 1 0 0],[1 0 0 0 0 1]},6};

test_fun = @(x,y) y + x;                element6{length(element6)+1}       = {test_fun,{b,b}, test_fun(b_val,b_val)     ,{[2 0 0 0 0 0],[0 2 0 0 0 0],[0 0 2 0 0 0],[0 0 0 2 0 0],[0 0 0 0 2 0],[0 0 0 0 0 2]},6};
test_fun = @(x,y) y + x;                element6{length(element6)+1}       = {test_fun,{b',b'},test_fun(b_val',b_val')  ,{[2 0 0 0 0 0];[0 2 0 0 0 0];[0 0 2 0 0 0];[0 0 0 2 0 0];[0 0 0 0 2 0];[0 0 0 0 0 2]},6};

test_fun = @(x,y) x + y;                element6{length(element6)+1}       = {test_fun,{B,B},   test_fun(B_val,B_val)   ,{[2 0 0 0 0 0],[0 0 2 0 0 0],[0 0 0 0 2 0];[0 2 0 0 0 0],[0 0 0 2 0 0],[0 0 0 0 0 2]},6};

%% MATRIX ARITHMETIC OPERATION - SEVENTH TEST SERIES MINUS(-)
test_fun = @(x) x - [2 7 6];            element7{1}                        = {test_fun,{x}, test_fun(x_val)          ,{[1 0 0 0 0],[1 0 0 0 0],[1 0 0 0 0]},5};
test_fun = @(x) x - [2;7;6];            element7{length(element7)+1}       = {test_fun,{x}, test_fun(x_val)          ,{[1 0 0 0 0];[1 0 0 0 0];[1 0 0 0 0]},5};
test_fun = @(x) x - [2 7 6;8 6 3];      element7{length(element7)+1}       = {test_fun,{x}, test_fun(x_val)          ,{[1 0 0 0 0],[1 0 0 0 0],[1 0 0 0 0];[1 0 0 0 0],[1 0 0 0 0],[1 0 0 0 0]},5};

test_fun = @(x) x - 8;                  element7{length(element7)+1}       = {test_fun,{b}, test_fun(b_val)          ,{[1 0 0 0 0 0],[0 1 0 0 0 0],[0 0 1 0 0 0],[0 0 0 1 0 0],[0 0 0 0 1 0],[0 0 0 0 0 1]},6};
test_fun = @(x) x - [2 7 6 8 6 3];      element7{length(element7)+1}       = {test_fun,{b}, test_fun(b_val)          ,{[1 0 0 0 0 0],[0 1 0 0 0 0],[0 0 1 0 0 0],[0 0 0 1 0 0],[0 0 0 0 1 0],[0 0 0 0 0 1]},6};

test_fun = @(x) x - 8;                  element7{length(element7)+1}       = {test_fun,{b'}, test_fun(b_val')        ,{[1 0 0 0 0 0];[0 1 0 0 0 0];[0 0 1 0 0 0];[0 0 0 1 0 0];[0 0 0 0 1 0];[0 0 0 0 0 1]},6};
test_fun = @(x) x - [2 7 6 8 6 3]';     element7{length(element7)+1}       = {test_fun,{b'}, test_fun(b_val')        ,{[1 0 0 0 0 0];[0 1 0 0 0 0];[0 0 1 0 0 0];[0 0 0 1 0 0];[0 0 0 0 1 0];[0 0 0 0 0 1]},6};

test_fun = @(x) x - 6.41;               element7{length(element7)+1}       = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix((ones(size(B_val)))),n_in_B};
test_fun = @(x) x - [2 7 6;8 6 3];      element7{length(element7)+1}       = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix(ones(size(B_val))),n_in_B};

%

test_fun = @(x) [2 7 6] - x;            element7{length(element7)+1}       = {test_fun,{x}, test_fun(x_val)          ,{[-1 0 0 0 0],[-1 0 0 0 0],[-1 0 0 0 0]},5};
test_fun = @(x) [2;7;6] - x;            element7{length(element7)+1}       = {test_fun,{x}, test_fun(x_val)          ,{[-1 0 0 0 0];[-1 0 0 0 0];[-1 0 0 0 0]},5};
test_fun = @(x) [2 7 6;8 6 3] - x;      element7{length(element7)+1}       = {test_fun,{x}, test_fun(x_val)          ,{[-1 0 0 0 0],[-1 0 0 0 0],[-1 0 0 0 0];[-1 0 0 0 0],[-1 0 0 0 0],[-1 0 0 0 0]},5};

test_fun = @(x) 8 - x;                  element7{length(element7)+1}       = {test_fun,{b}, test_fun(b_val)          ,{[-1 0 0 0 0 0],[0 -1 0 0 0 0],[0 0 -1 0 0 0],[0 0 0 -1 0 0],[0 0 0 0 -1 0],[0 0 0 0 0 -1]},6};
test_fun = @(x) [2 7 6 8 6 3] - x;      element7{length(element7)+1}       = {test_fun,{b}, test_fun(b_val)          ,{[-1 0 0 0 0 0],[0 -1 0 0 0 0],[0 0 -1 0 0 0],[0 0 0 -1 0 0],[0 0 0 0 -1 0],[0 0 0 0 0 -1]},6};

test_fun = @(x) 8 - x;                  element7{length(element7)+1}       = {test_fun,{b'}, test_fun(b_val')        ,{[-1 0 0 0 0 0];[0 -1 0 0 0 0];[0 0 -1 0 0 0];[0 0 0 -1 0 0];[0 0 0 0 -1 0];[0 0 0 0 0 -1]},6};
test_fun = @(x) [2 7 6 8 6 3]' - x;     element7{length(element7)+1}       = {test_fun,{b'}, test_fun(b_val')        ,{[-1 0 0 0 0 0];[0 -1 0 0 0 0];[0 0 -1 0 0 0];[0 0 0 -1 0 0];[0 0 0 0 -1 0];[0 0 0 0 0 -1]},6};

test_fun = @(x) 8 - x;                  element7{length(element7)+1}       = {test_fun,{b'}, test_fun(b_val')        ,{[-1 0 0 0 0 0];[0 -1 0 0 0 0];[0 0 -1 0 0 0];[0 0 0 -1 0 0];[0 0 0 0 -1 0];[0 0 0 0 0 -1]},6};
test_fun = @(x) [2 7 6 8 6 3]' - x;     element7{length(element7)+1}       = {test_fun,{b'}, test_fun(b_val')        ,{[-1 0 0 0 0 0];[0 -1 0 0 0 0];[0 0 -1 0 0 0];[0 0 0 -1 0 0];[0 0 0 0 -1 0];[0 0 0 0 0 -1]},6};

test_fun = @(x) 6.41 - x;               element7{length(element7)+1}       = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix((-ones(size(B_val)))),n_in_B};
test_fun = @(x) [2 7 6;8 6 3] - x;      element7{length(element7)+1}       = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix(-ones(size(B_val))),n_in_B};

%

test_fun = @(x,y) x - y;                element7{length(element7)+1}       = {test_fun,{x,y},   test_fun(x_val,y_val)   ,{[1 -1 0 0 0]},5};
test_fun = @(x,y) x - y;                element7{length(element7)+1}       = {test_fun,{x,yy},  test_fun(x_val,yy_val)  ,{[1 -1 0 0 0],[1 0 -1 0 0]},5};
test_fun = @(x,y) x - y;                element7{length(element7)+1}       = {test_fun,{x,yy'}, test_fun(x_val,yy_val') ,{[1 -1 0 0 0];[1 0 -1 0 0]},5};
test_fun = @(x,y) x - y;                element7{length(element7)+1}       = {test_fun,{x,B},   test_fun(x_val,B_val)   ,{[0 0 0 0 0 0],[1 0 -1 0 0 0],[1 0 0 0 -1 0];[1 -1 0 0 0 0],[1 0 0 -1 0 0],[1 0 0 0 0 -1]},6};

test_fun = @(x,y) x - y;                element7{length(element7)+1}       = {test_fun,{yy,x},  test_fun(yy_val,x_val)  ,{[-1 1 0 0 0],[-1 0 1 0 0]},5};
test_fun = @(x,y) x - y;                element7{length(element7)+1}       = {test_fun,{yy',x}, test_fun(yy_val',x_val) ,{[-1 1 0 0 0];[-1 0 1 0 0]},5};
test_fun = @(x,y) x - y;                element7{length(element7)+1}       = {test_fun,{B,x},   test_fun(B_val,x_val)   ,{[0 0 0 0 0 0],[-1 0 1 0 0 0],[-1 0 0 0 1 0];[-1 1 0 0 0 0],[-1 0 0 1 0 0],[-1 0 0 0 0 1]},6};

test_fun = @(x,y) y - x;                element7{length(element7)+1}       = {test_fun,{b,b}, test_fun(b_val,b_val)     ,{[0 0 0 0 0 0],[0 0 0 0 0 0],[0 0 0 0 0 0],[0 0 0 0 0 0],[0 0 0 0 0 0],[0 0 0 0 0 0]},6};
test_fun = @(x,y) y - x;                element7{length(element7)+1}       = {test_fun,{b',b'},test_fun(b_val',b_val')  ,{[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0]},6};

test_fun = @(x,y) x - y;                element6{length(element6)+1}       = {test_fun,{B,B},   test_fun(B_val,B_val)   ,{[0 0 0 0 0 0],[0 0 0 0 0 0],[0 0 0 0 0 0];[0 0 0 0 0 0],[0 0 0 0 0 0],[0 0 0 0 0 0]},6};

%% MATRIX ARITHMETIC OPERATION - EIGHT TEST SERIES TIMES(.*)
test_fun = @(x) x .* [2 7 6];            element8{1}                       = {test_fun,{x}, test_fun(x_val)          ,{[2 0 0 0 0],[7 0 0 0 0],[6 0 0 0 0]},5};
test_fun = @(x) x .* [2;7;6];            element8{length(element8)+1}      = {test_fun,{x}, test_fun(x_val)          ,{[2 0 0 0 0];[7 0 0 0 0];[6 0 0 0 0]},5};
test_fun = @(x) x .* [2 7 6;8 6 3];      element8{length(element8)+1}      = {test_fun,{x}, test_fun(x_val)          ,{[2 0 0 0 0],[7 0 0 0 0],[6 0 0 0 0];[8 0 0 0 0],[6 0 0 0 0],[3 0 0 0 0]},5};

test_fun = @(x) x .* 8;                  element8{length(element8)+1}      = {test_fun,{b}, test_fun(b_val)          ,{[8 0 0 0 0 0],[0 8 0 0 0 0],[0 0 8 0 0 0],[0 0 0 8 0 0],[0 0 0 0 8 0],[0 0 0 0 0 8]},6};
test_fun = @(x) x .* [2 7 6 8 6 3];      element8{length(element8)+1}      = {test_fun,{b}, test_fun(b_val)          ,{[2 0 0 0 0 0],[0 7 0 0 0 0],[0 0 6 0 0 0],[0 0 0 8 0 0],[0 0 0 0 6 0],[0 0 0 0 0 3]},6};

test_fun = @(x) x .* 8;                  element8{length(element8)+1}      = {test_fun,{b'}, test_fun(b_val')        ,{[8 0 0 0 0 0];[0 8 0 0 0 0];[0 0 8 0 0 0];[0 0 0 8 0 0];[0 0 0 0 8 0];[0 0 0 0 0 8]},6};
test_fun = @(x) x .* [2 7 6 8 6 3]';     element8{length(element8)+1}      = {test_fun,{b'}, test_fun(b_val')        ,{[2 0 0 0 0 0];[0 7 0 0 0 0];[0 0 6 0 0 0];[0 0 0 8 0 0];[0 0 0 0 6 0];[0 0 0 0 0 3]},6};

test_fun = @(x) x .* 6.41;               element8{length(element8)+1}      = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix(6.41*(ones(size(B_val)))),n_in_B};
test_fun = @(x) x .* [2 7 6;8 6 3];      element8{length(element8)+1}      = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix([2 7 6;8 6 3]),n_in_B};

%
test_fun = @(x) [2 7 6] .* x;            element8{length(element8)+1}      = {test_fun,{x}, test_fun(x_val)          ,{[2 0 0 0 0],[7 0 0 0 0],[6 0 0 0 0]},5};
test_fun = @(x) [2;7;6] .* x;            element8{length(element8)+1}      = {test_fun,{x}, test_fun(x_val)          ,{[2 0 0 0 0];[7 0 0 0 0];[6 0 0 0 0]},5};
test_fun = @(x) [2 7 6;8 6 3] .* x;      element8{length(element8)+1}      = {test_fun,{x}, test_fun(x_val)          ,{[2 0 0 0 0],[7 0 0 0 0],[6 0 0 0 0];[8 0 0 0 0],[6 0 0 0 0],[3 0 0 0 0]},5};

test_fun = @(x) 8 .* x;                  element8{length(element8)+1}      = {test_fun,{b}, test_fun(b_val)          ,{[8 0 0 0 0 0],[0 8 0 0 0 0],[0 0 8 0 0 0],[0 0 0 8 0 0],[0 0 0 0 8 0],[0 0 0 0 0 8]},6};
test_fun = @(x) [2 7 6 8 6 3] .* x;      element8{length(element8)+1}      = {test_fun,{b}, test_fun(b_val)          ,{[2 0 0 0 0 0],[0 7 0 0 0 0],[0 0 6 0 0 0],[0 0 0 8 0 0],[0 0 0 0 6 0],[0 0 0 0 0 3]},6};

test_fun = @(x) 8 .* x;                  element8{length(element8)+1}      = {test_fun,{b'}, test_fun(b_val')        ,{[8 0 0 0 0 0];[0 8 0 0 0 0];[0 0 8 0 0 0];[0 0 0 8 0 0];[0 0 0 0 8 0];[0 0 0 0 0 8]},6};
test_fun = @(x) [2 7 6 8 6 3]' .* x;     element8{length(element8)+1}      = {test_fun,{b'}, test_fun(b_val')        ,{[2 0 0 0 0 0];[0 7 0 0 0 0];[0 0 6 0 0 0];[0 0 0 8 0 0];[0 0 0 0 6 0];[0 0 0 0 0 3]},6};

test_fun = @(x) 6.41 .* x;               element8{length(element8)+1}      = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix(6.41*(ones(size(B_val)))),n_in_B};
test_fun = @(x) [2 7 6;8 6 3] .* x;      element8{length(element8)+1}      = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix([2 7 6;8 6 3]),n_in_B};

%18up 19 down

test_fun = @(x,y) x .* y;                element8{length(element8)+1}      = {test_fun,{x,y},   test_fun(x_val,y_val)   ,{[y_val x_val 0 0 0]},5};
test_fun = @(x,y) x .* y;                element8{length(element8)+1}      = {test_fun,{x,yy},  test_fun(x_val,yy_val)  ,{[yy_val(1) x_val 0 0 0],[yy_val(2) 0 x_val 0 0]},5};
test_fun = @(x,y) x .* y;                element8{length(element8)+1}      = {test_fun,{x,yy'}, test_fun(x_val,yy_val') ,{[yy_val(1) x_val 0 0 0];[yy_val(2) 0 x_val 0 0]},5};
test_fun = @(x,y) x .* y;                element8{length(element8)+1}      = {test_fun,{x,B},   test_fun(x_val,B_val)   ,{[x_val+B_val(1,1) 0 0 0 0 0],[B_val(1,2) 0 x_val 0 0 0],[B_val(1,3) 0 0 0 x_val 0];[B_val(2,1) x_val 0 0 0 0],[B_val(2,2) 0 0 x_val 0 0],[B_val(2,3) 0 0 0 0 x_val]},6};

test_fun = @(x,y) x .* y;                element8{length(element8)+1}      = {test_fun,{yy,x},  test_fun(yy_val,x_val)  ,{[yy_val(1) x_val 0 0 0],[yy_val(2) 0 x_val 0 0]},5};
test_fun = @(x,y) x .* y;                element8{length(element8)+1}      = {test_fun,{yy',x}, test_fun(yy_val',x_val) ,{[yy_val(1) x_val 0 0 0];[yy_val(2) 0 x_val 0 0]},5};
test_fun = @(x,y) x .* y;                element8{length(element8)+1}      = {test_fun,{B,x},   test_fun(B_val,x_val)   ,{[x_val+B_val(1,1) 0 0 0 0 0],[B_val(1,2) 0 x_val 0 0 0],[B_val(1,3) 0 0 0 x_val 0];[B_val(2,1) x_val 0 0 0 0],[B_val(2,2) 0 0 x_val 0 0],[B_val(2,3) 0 0 0 0 x_val]},6};

test_fun = @(x,y) y .* x;                element8{length(element8)+1}      = {test_fun,{b,b}, test_fun(b_val,b_val)     ,{[2.*b_val(1) 0 0 0 0 0],[0 2.*b_val(2) 0 0 0 0],[0 0 2.*b_val(3) 0 0 0],[0 0 0 2.*b_val(4) 0 0],[0 0 0 0 2.*b_val(5) 0],[0 0 0 0 0 2.*b_val(6)]},6};
test_fun = @(x,y) y .* x;                element8{length(element8)+1}      = {test_fun,{b',b'},test_fun(b_val',b_val')  ,{[2.*b_val(1) 0 0 0 0 0];[0 2.*b_val(2) 0 0 0 0];[0 0 2.*b_val(3) 0 0 0];[0 0 0 2.*b_val(4) 0 0];[0 0 0 0 2.*b_val(5) 0];[0 0 0 0 0 2.*b_val(6)]},6};

test_fun = @(x,y) x .* y;                element8{length(element8)+1}      = {test_fun,{B,B},   test_fun(B_val,B_val)   ,{[2.*B_val(1,1) 0 0 0 0 0],[0 0 2.*B_val(1,2) 0 0 0],[0 0 0 0 2.*B_val(1,3) 0];[0 2.*B_val(2,1) 0 0 0 0],[0 0 0 2.*B_val(2,2) 0 0],[0 0 0 0 0 2.*B_val(2,3)]},6};

%% MATRIX ARITHMETIC OPERATION - NINTH TEST SERIES POWER(.^)
test_fun = @(x) x .^ 6;                  element9{1}                       = {test_fun,{x}, test_fun(x_val)          ,{[6.*x_val.^5 0 0 0 0]},5};
test_fun = @(x) x .^ [2 7 6];            element9{length(element9)+1}      = {test_fun,{x}, test_fun(x_val)          ,{[2.*x_val 0 0 0 0],[7.*x_val.^6 0 0 0 0],[6.*x_val.^5 0 0 0 0]},5};
test_fun = @(x) x .^ [2;7;6];            element9{length(element9)+1}      = {test_fun,{x}, test_fun(x_val)          ,{[2.*x_val 0 0 0 0];[7.*x_val.^6 0 0 0 0];[6.*x_val.^5 0 0 0 0]},5};
test_fun = @(x) x .^ [2 7 6;8 6 3];      element9{length(element9)+1}      = {test_fun,{x}, test_fun(x_val)          ,{[2.*x_val 0 0 0 0],[7.*x_val.^6 0 0 0 0],[6.*x_val.^5 0 0 0 0];[8.*x_val.^7 0 0 0 0],[6.*x_val.^5 0 0 0 0],[3.*x_val.^2 0 0 0 0]},5};

test_fun = @(x) x .^ 8;                  element9{length(element9)+1}      = {test_fun,{b}, test_fun(b_val)          ,{[8.*b_val(1).^7 0 0 0 0 0],[0 8.*b_val(2).^7 0 0 0 0],[0 0 8.*b_val(3).^7 0 0 0],[0 0 0 8.*b_val(4).^7 0 0],[0 0 0 0 8.*b_val(5).^7 0],[0 0 0 0 0 8.*b_val(6).^7]},6};
test_fun = @(x) x .^ [2 7 6 8 6 3];      element9{length(element9)+1}      = {test_fun,{b}, test_fun(b_val)          ,{[2.*b_val(1).^1 0 0 0 0 0],[0 7.*b_val(2).^6 0 0 0 0],[0 0 6.*b_val(3).^5 0 0 0],[0 0 0 8.*b_val(4).^7 0 0],[0 0 0 0 6.*b_val(5).^5 0],[0 0 0 0 0 3.*b_val(6).^2]},6};

test_fun = @(x) x .^ 8;                  element9{length(element9)+1}      = {test_fun,{b'}, test_fun(b_val')        ,{[8.*b_val(1).^7 0 0 0 0 0];[0 8.*b_val(2).^7 0 0 0 0];[0 0 8.*b_val(3).^7 0 0 0];[0 0 0 8.*b_val(4).^7 0 0];[0 0 0 0 8.*b_val(5).^7 0];[0 0 0 0 0 8.*b_val(6).^7]},6};
test_fun = @(x) x .^ ([2 7 6 8 6 3]');   element9{length(element9)+1}      = {test_fun,{b'}, test_fun(b_val')        ,{[2.*b_val(1).^1 0 0 0 0 0];[0 7.*b_val(2).^6 0 0 0 0];[0 0 6.*b_val(3).^5 0 0 0];[0 0 0 8.*b_val(4).^7 0 0];[0 0 0 0 6.*b_val(5).^5 0];[0 0 0 0 0 3.*b_val(6).^2]},6};

test_fun = @(x) x .^ 6.41;               element9{length(element9)+1}      = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix(6.41*(B_val.^5.41)),n_in_B};
test_fun = @(x) x .^ [2 7 6;8 6 3];      element9{length(element9)+1}      = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix([2 7 6;8 6 3].*B_val.^([2 7 6;8 6 3]-1)),n_in_B};

%

test_fun = @(x) 6 .^ x;                  element9{length(element9)+1}      = {test_fun,{x}, test_fun(x_val)          ,{[log(6).*6.^x_val 0 0 0 0]},5};
test_fun = @(x) [2 7 6] .^ x;            element9{length(element9)+1}      = {test_fun,{x}, test_fun(x_val)          ,{[log(2).*2.^x_val 0 0 0 0],[log(7).*7.^x_val 0 0 0 0],[log(6).*6.^x_val 0 0 0 0]},5};
test_fun = @(x) [2;7;6] .^ x;            element9{length(element9)+1}      = {test_fun,{x}, test_fun(x_val)          ,{[log(2).*2.^x_val 0 0 0 0];[log(7).*7.^x_val 0 0 0 0];[log(6).*6.^x_val 0 0 0 0]},5};
test_fun = @(x) [2 7 6;8 6 3] .^ x;      element9{length(element9)+1}      = {test_fun,{x}, test_fun(x_val)          ,{[log(2).*2.^x_val 0 0 0 0],[log(7).*7.^x_val 0 0 0 0],[log(6).*6.^x_val 0 0 0 0];[log(8).*8.^x_val 0 0 0 0],[log(6).*6.^x_val 0 0 0 0],[log(3).*3.^x_val 0 0 0 0]},5};

test_fun = @(x) 8 .^ x;                  element9{length(element9)+1}      = {test_fun,{b}, test_fun(b_val)          ,{[log(8).*8.^b_val(1) 0 0 0 0 0],[0 log(8).*8.^b_val(2) 0 0 0 0],[0 0 log(8).*8.^b_val(3) 0 0 0],[0 0 0 log(8).*8.^b_val(4) 0 0],[0 0 0 0 log(8).*8.^b_val(5) 0],[0 0 0 0 0 log(8).*8.^b_val(6)]},6};
test_fun = @(x) [2 7 6 8 6 3] .^ x;      element9{length(element9)+1}      = {test_fun,{b}, test_fun(b_val)          ,{[log(2).*2.^b_val(1) 0 0 0 0 0],[0 log(7).*7.^b_val(2) 0 0 0 0],[0 0 log(6).*6.^b_val(3) 0 0 0],[0 0 0 log(8).*8.^b_val(4) 0 0],[0 0 0 0 log(6).*6.^b_val(5) 0],[0 0 0 0 0 log(3).*3.^b_val(6)]},6};

test_fun = @(x) 8 .^ x;                  element9{length(element9)+1}      = {test_fun,{b'}, test_fun(b_val')        ,{[log(8).*8.^b_val(1) 0 0 0 0 0];[0 log(8).*8.^b_val(2) 0 0 0 0];[0 0 log(8).*8.^b_val(3) 0 0 0];[0 0 0 log(8).*8.^b_val(4) 0 0];[0 0 0 0 log(8).*8.^b_val(5) 0];[0 0 0 0 0 log(8).*8.^b_val(6)]},6};
test_fun = @(x) ([2 7 6 8 6 3]') .^ x;   element9{length(element9)+1}      = {test_fun,{b'}, test_fun(b_val')        ,{[log(2).*2.^b_val(1) 0 0 0 0 0];[0 log(7).*7.^b_val(2) 0 0 0 0];[0 0 log(6).*6.^b_val(3) 0 0 0];[0 0 0 log(8).*8.^b_val(4) 0 0];[0 0 0 0 log(6).*6.^b_val(5) 0];[0 0 0 0 0 log(3).*3.^b_val(6)]},6};

test_fun = @(x) 6.41 .^ x;               element9{length(element9)+1}      = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix(log(6.41)*6.41.^B_val),n_in_B};
test_fun = @(x) [2 7 6;8 6 3] .^ x;      element9{length(element9)+1}      = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix(log([2 7 6;8 6 3]).*[2 7 6;8 6 3].^B_val),n_in_B};
 
%%%20 up and 21 down

test_fun = @(x,y) x .^ y;                element9{length(element9)+1}      = {test_fun,{x,y},   test_fun(x_val,y_val)   ,{[y_val.*x_val.^(y_val-1) log(x_val).*x_val.^y_val 0 0 0]},5};
test_fun = @(x,y) x .^ y;                element9{length(element9)+1}      = {test_fun,{x,yy},  test_fun(x_val,yy_val)  ,{[yy_val(1).*x_val.^(yy_val(1)-1) log(x_val).*x_val.^yy_val(1) 0 0 0],[yy_val(2).*x_val.^(yy_val(2)-1) 0 log(x_val).*x_val.^yy_val(2) 0 0]},5};
test_fun = @(x,y) x .^ y;                element9{length(element9)+1}      = {test_fun,{x,yy'}, test_fun(x_val,yy_val') ,{[yy_val(1).*x_val.^(yy_val(1)-1) log(x_val).*x_val.^yy_val(1) 0 0 0];[yy_val(2).*x_val.^(yy_val(2)-1) 0 log(x_val).*x_val.^yy_val(2) 0 0]},5};
test_fun = @(x,y) x .^ y;                element9{length(element9)+1}      = {test_fun,{x,B},   test_fun(x_val,B_val)   ,{[B_val(1,1).*x_val.^(B_val(1,1)-1)+log(x_val).*x_val.^B_val(1,1) 0 0 0 0 0],[B_val(1,2).*x_val.^(B_val(1,2)-1) 0 log(x_val).*x_val.^B_val(1,2) 0 0 0],[B_val(1,3).*x_val.^(B_val(1,3)-1) 0 0 0 log(x_val).*x_val.^B_val(1,3) 0];[B_val(2,1).*x_val.^(B_val(2,1)-1) log(x_val).*x_val.^B_val(2,1) 0 0 0 0],[B_val(2,2).*x_val.^(B_val(2,2)-1) 0 0 log(x_val).*x_val.^B_val(2,2) 0 0],[B_val(2,3).*x_val.^(B_val(2,3)-1) 0 0 0 0 log(x_val).*x_val.^B_val(2,3)]},6};

test_fun = @(x,y) x .^ y;                element9{length(element9)+1}      = {test_fun,{yy,x},  test_fun(yy_val,x_val)  ,{[log(yy_val(1)).*yy_val(1).^x_val x_val.*yy_val(1).^(x_val-1) 0 0 0],[log(yy_val(2)).*yy_val(2).^x_val 0 x_val.*yy_val(2).^(x_val-1) 0 0]},5};
test_fun = @(x,y) x .^ y;                element9{length(element9)+1}      = {test_fun,{yy',x}, test_fun(yy_val',x_val) ,{[log(yy_val(1)).*yy_val(1).^x_val x_val.*yy_val(1).^(x_val-1) 0 0 0];[log(yy_val(2)).*yy_val(2).^x_val 0 x_val.*yy_val(2).^(x_val-1) 0 0]},5};
test_fun = @(x,y) x .^ y;                element9{length(element9)+1}      = {test_fun,{B,x},   test_fun(B_val,x_val)   ,{[log(B_val(1,1)).*B_val(1,1).^x_val + x_val.*B_val(1,1).^(x_val-1) 0 0 0 0 0],[log(B_val(1,2)).*B_val(1,2).^x_val 0 x_val.*B_val(1,2).^(x_val-1) 0 0 0],[log(B_val(1,3)).*B_val(1,3).^x_val 0 0 0 x_val.*B_val(1,3).^(x_val-1) 0];[log(B_val(2,1)).*B_val(2,1).^x_val x_val.*B_val(2,1).^(x_val-1) 0 0 0 0],[log(B_val(2,2)).*B_val(2,2).^x_val 0 0 x_val.*B_val(2,2).^(x_val-1) 0 0],[log(B_val(2,3)).*B_val(2,3).^x_val 0 0 0 0 x_val.*B_val(2,3).^(x_val-1)]},6};
% 
test_fun = @(x,y) y .^ x;                element9{length(element9)+1}      = {test_fun,{b,b}, test_fun(b_val,b_val)     ,{[log(b_val(1)).*b_val(1).^b_val(1) + b_val(1).*b_val(1).^(b_val(1)-1) 0 0 0 0 0],[0 log(b_val(2)).*b_val(2).^b_val(2) + b_val(2).*b_val(2).^(b_val(2)-1) 0 0 0 0],[0 0 log(b_val(3)).*b_val(3).^b_val(3) + b_val(3).*b_val(3).^(b_val(3)-1) 0 0 0],[0 0 0 log(b_val(4)).*b_val(4).^b_val(4) + b_val(4).*b_val(4).^(b_val(4)-1) 0 0],[0 0 0 0 log(b_val(5)).*b_val(5).^b_val(5) + b_val(5).*b_val(5).^(b_val(5)-1) 0],[0 0 0 0 0 log(b_val(6)).*b_val(6).^b_val(6) + b_val(6).*b_val(6).^(b_val(6)-1)]},6};
test_fun = @(x,y) y .^ x;                element9{length(element9)+1}      = {test_fun,{b',b'},test_fun(b_val',b_val')  ,{[log(b_val(1)).*b_val(1).^b_val(1) + b_val(1).*b_val(1).^(b_val(1)-1) 0 0 0 0 0];[0 log(b_val(2)).*b_val(2).^b_val(2) + b_val(2).*b_val(2).^(b_val(2)-1) 0 0 0 0];[0 0 log(b_val(3)).*b_val(3).^b_val(3) + b_val(3).*b_val(3).^(b_val(3)-1) 0 0 0];[0 0 0 log(b_val(4)).*b_val(4).^b_val(4) + b_val(4).*b_val(4).^(b_val(4)-1) 0 0];[0 0 0 0 log(b_val(5)).*b_val(5).^b_val(5) + b_val(5).*b_val(5).^(b_val(5)-1) 0];[0 0 0 0 0 log(b_val(6)).*b_val(6).^b_val(6) + b_val(6).*b_val(6).^(b_val(6)-1)]},6};
% 
test_fun = @(x,y) x .^ y;                element9{length(element9)+1}      = {test_fun,{B,B},   test_fun(B_val,B_val)   ,{[log(B_val(1,1)).*B_val(1,1).^B_val(1,1) + B_val(1,1).*B_val(1,1).^(B_val(1,1)-1) 0 0 0 0 0],[0 0 log(B_val(1,2)).*B_val(1,2).^B_val(1,2) + B_val(1,2).*B_val(1,2).^(B_val(1,2)-1) 0 0 0],[0 0 0 0 log(B_val(1,3)).*B_val(1,3).^B_val(1,3) + B_val(1,3).*B_val(1,3).^(B_val(1,3)-1) 0];[0 log(B_val(2,1)).*B_val(2,1).^B_val(2,1) + B_val(2,1).*B_val(2,1).^(B_val(2,1)-1) 0 0 0 0],[0 0 0 log(B_val(2,2)).*B_val(2,2).^B_val(2,2) + B_val(2,2).*B_val(2,2).^(B_val(2,2)-1) 0 0],[0 0 0 0 0 log(B_val(2,3)).*B_val(2,3).^B_val(2,3) + B_val(2,3).*B_val(2,3).^(B_val(2,3)-1)]},6};

%% MATRIX ARITHMETIC OPERATION - TENTH TEST SERIES DIVIDE(./)
test_fun = @(x) x ./ [2 7 6];            element10{1}                      = {test_fun,{x}, test_fun(x_val)          ,{[1/2 0 0 0 0],[1/7 0 0 0 0],[1/6 0 0 0 0]},5};
test_fun = @(x) x ./ [2;7;6];            element10{length(element10)+1}    = {test_fun,{x}, test_fun(x_val)          ,{[1/2 0 0 0 0];[1/7 0 0 0 0];[1/6 0 0 0 0]},5};
test_fun = @(x) x ./ [2 7 6;8 6 3];      element10{length(element10)+1}    = {test_fun,{x}, test_fun(x_val)          ,{[1/2 0 0 0 0],[1/7 0 0 0 0],[1/6 0 0 0 0];[1/8 0 0 0 0],[1/6 0 0 0 0],[1/3 0 0 0 0]},5};

test_fun = @(x) x ./ 8;                  element10{length(element10)+1}    = {test_fun,{b}, test_fun(b_val)          ,{[1/8 0 0 0 0 0],[0 1/8 0 0 0 0],[0 0 1/8 0 0 0],[0 0 0 1/8 0 0],[0 0 0 0 1/8 0],[0 0 0 0 0 1/8]},6};
test_fun = @(x) x ./ [2 7 6 8 6 3];      element10{length(element10)+1}    = {test_fun,{b}, test_fun(b_val)          ,{[1/2 0 0 0 0 0],[0 1/7 0 0 0 0],[0 0 1/6 0 0 0],[0 0 0 1/8 0 0],[0 0 0 0 1/6 0],[0 0 0 0 0 1/3]},6};

test_fun = @(x) x ./ 8;                  element10{length(element10)+1}    = {test_fun,{b'}, test_fun(b_val')        ,{[1/8 0 0 0 0 0];[0 1/8 0 0 0 0];[0 0 1/8 0 0 0];[0 0 0 1/8 0 0];[0 0 0 0 1/8 0];[0 0 0 0 0 1/8]},6};
test_fun = @(x) x ./ [2 7 6 8 6 3]';     element10{length(element10)+1}    = {test_fun,{b'}, test_fun(b_val')        ,{[1/2 0 0 0 0 0];[0 1/7 0 0 0 0];[0 0 1/6 0 0 0];[0 0 0 1/8 0 0];[0 0 0 0 1/6 0];[0 0 0 0 0 1/3]},6};

test_fun = @(x) x ./ 6.41;               element10{length(element10)+1}    = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix(1/6.41*(ones(size(B_val)))),n_in_B};
test_fun = @(x) x ./ [2 7 6;8 6 3];      element10{length(element10)+1}    = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix([1/2 1/7 1/6;1/8 1/6 1/3]),n_in_B};

%
test_fun = @(x) [2 7 6] ./ x;            element10{length(element10)+1}    = {test_fun,{x}, test_fun(x_val)          ,{[-2/(x_val^2) 0 0 0 0],[-7/(x_val^2) 0 0 0 0],[-6/(x_val^2) 0 0 0 0],[-7/(x_val^2) 0 0 0 0],[-6/(x_val^2) 0 0 0 0]},5};
test_fun = @(x) [2;7;6] ./ x;            element10{length(element10)+1}    = {test_fun,{x}, test_fun(x_val)          ,{[-2/(x_val^2) 0 0 0 0];[-7/(x_val^2) 0 0 0 0];[-6/(x_val^2) 0 0 0 0]},5};
test_fun = @(x) [2 7 6;8 6 3] ./ x;      element10{length(element10)+1}    = {test_fun,{x}, test_fun(x_val)          ,{[-2/(x_val^2) 0 0 0 0],[-7/(x_val^2) 0 0 0 0],[-6/(x_val^2) 0 0 0 0];[-8/(x_val^2) 0 0 0 0],[-6/(x_val^2) 0 0 0 0],[-3/(x_val^2) 0 0 0 0]},5};

test_fun = @(x) 8 ./ x;                  element10{length(element10)+1}    = {test_fun,{b}, test_fun(b_val)          ,{[-8/(b_val(1)^2) 0 0 0 0 0],[0 -8/(b_val(2)^2) 0 0 0 0],[0 0 -8/(b_val(3)^2) 0 0 0],[0 0 0 -8/(b_val(4)^2) 0 0],[0 0 0 0 -8/(b_val(5)^2) 0],[0 0 0 0 0 -8/(b_val(6)^2)]},6};
test_fun = @(x) [2 7 6 8 6 3] ./ x;      element10{length(element10)+1}    = {test_fun,{b}, test_fun(b_val)          ,{[-2/(b_val(1)^2) 0 0 0 0 0],[0 -7/(b_val(2)^2) 0 0 0 0],[0 0 -6/(b_val(3)^2) 0 0 0],[0 0 0 -8/(b_val(4)^2) 0 0],[0 0 0 0 -6/(b_val(5)^2) 0],[0 0 0 0 0 -3/(b_val(6)^2)]},6};

test_fun = @(x) 8 ./ x;                  element10{length(element10)+1}    = {test_fun,{b'}, test_fun(b_val')        ,{[-8/(b_val(1)^2) 0 0 0 0 0];[0 -8/(b_val(2)^2) 0 0 0 0];[0 0 -8/(b_val(3)^2) 0 0 0];[0 0 0 -8/(b_val(4)^2) 0 0];[0 0 0 0 -8/(b_val(5)^2) 0];[0 0 0 0 0 -8/(b_val(6)^2)]},6};
test_fun = @(x) [2 7 6 8 6 3]' ./ x;     element10{length(element10)+1}    = {test_fun,{b'}, test_fun(b_val')        ,{[-2/(b_val(1)^2) 0 0 0 0 0];[0 -7/(b_val(2)^2) 0 0 0 0];[0 0 -6/(b_val(3)^2) 0 0 0];[0 0 0 -8/(b_val(4)^2) 0 0];[0 0 0 0 -6/(b_val(5)^2) 0];[0 0 0 0 0 -3/(b_val(6)^2)]},6};

test_fun = @(x) 6.41 ./ x;               element10{length(element10)+1}    = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix(-6.41./(B_val.^2)),n_in_B};
test_fun = @(x) [2 7 6;8 6 3] ./ x;      element10{length(element10)+1}    = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix(-[2 7 6;8 6 3]./(B_val.^2)),n_in_B};

%18up 19 down

test_fun = @(x,y) x ./ y;                element10{length(element10)+1}    = {test_fun,{x,y},   test_fun(x_val,y_val)   ,{[1/y_val -x_val/(y_val^2) 0 0 0]},5};
test_fun = @(x,y) x ./ y;                element10{length(element10)+1}    = {test_fun,{x,yy},  test_fun(x_val,yy_val)  ,{[1/yy_val(1) -x_val/(yy_val(1)^2) 0 0 0],[1/yy_val(2) 0 -x_val/(yy_val(2)^2) 0 0]},5};
test_fun = @(x,y) x ./ y;                element10{length(element10)+1}    = {test_fun,{x,yy'}, test_fun(x_val,yy_val') ,{[1/yy_val(1) -x_val/(yy_val(1)^2) 0 0 0];[1/yy_val(2) 0 -x_val/(yy_val(2)^2) 0 0]},5};
test_fun = @(x,y) x ./ y;                element10{length(element10)+1}    = {test_fun,{x,B},   test_fun(x_val,B_val)   ,{[1/B_val(1,1)-x_val/(B_val(1,1)^2) 0 0 0 0 0],[1/B_val(1,2) 0 -x_val/(B_val(1,2)^2) 0 0 0],[1/B_val(1,3) 0 0 0 -x_val/(B_val(1,3)^2) 0];[1/B_val(2,1) -x_val/(B_val(2,1)^2) 0 0 0 0],[1/B_val(2,2) 0 0 -x_val/(B_val(2,2)^2) 0 0],[1/B_val(2,3) 0 0 0 0 -x_val/(B_val(2,3)^2)]},6};

test_fun = @(x,y) x ./ y;                element10{length(element10)+1}    = {test_fun,{yy,x},  test_fun(yy_val,x_val)  ,{[-yy_val(1)/(x_val^2) 1/x_val 0 0 0],[-yy_val(2)/(x_val^2) 0 1/x_val 0 0]},5};
test_fun = @(x,y) x ./ y;                element10{length(element10)+1}    = {test_fun,{yy',x}, test_fun(yy_val',x_val) ,{[-yy_val(1)/(x_val^2) 1/x_val 0 0 0];[-yy_val(2)/(x_val^2) 0 1/x_val 0 0]},5};
test_fun = @(x,y) x ./ y;                element10{length(element10)+1}    = {test_fun,{B,x},   test_fun(B_val,x_val)   ,{[1/x_val-B_val(1,1)/(x_val^2) 0 0 0 0 0],[-B_val(1,2)/(x_val^2) 0 1/x_val 0 0 0],[-B_val(1,3)/(x_val^2) 0 0 0 1/x_val 0];[-B_val(2,1)/(x_val^2) 1/x_val 0 0 0 0],[-B_val(2,2)/(x_val^2) 0 0 1/x_val 0 0],[-B_val(2,3)/(x_val^2) 0 0 0 0 1/x_val]},6};

test_fun = @(x,y) y ./ x;                element10{length(element10)+1}    = {test_fun,{b,b}, test_fun(b_val,b_val)     ,{[0 0 0 0 0 0],[0 0 0 0 0 0],[0 0 0 0 0 0],[0 0 0 0 0 0],[0 0 0 0 0 0],[0 0 0 0 0 0]},6};
test_fun = @(x,y) y ./ x;                element10{length(element10)+1}    = {test_fun,{b',b'},test_fun(b_val',b_val')  ,{[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0]},6};

test_fun = @(x,y) x ./ y;                element10{length(element10)+1}    = {test_fun,{B,B},   test_fun(B_val,B_val)   ,{[0 0 0 0 0 0],[0 0 0 0 0 0],[0 0 0 0 0 0];[0 0 0 0 0 0],[0 0 0 0 0 0],[0 0 0 0 0 0]},6};

%% MATRIX ARITHMETIC OPERATION - ELEVENTH TEST SERIES RDIVIDE(/)

test_fun = @(x) x / 8;                   element11{1}                      = {test_fun,{b}, test_fun(b_val)             ,{[1/8 0 0 0 0 0],[0 1/8 0 0 0 0],[0 0 1/8 0 0 0],[0 0 0 1/8 0 0],[0 0 0 0 1/8 0],[0 0 0 0 0 1/8]},6};
test_fun = @(x) x / 8;                   element11{length(element11)+1}    = {test_fun,{b'}, test_fun(b_val')           ,{[1/8 0 0 0 0 0];[0 1/8 0 0 0 0];[0 0 1/8 0 0 0];[0 0 0 1/8 0 0];[0 0 0 0 1/8 0];[0 0 0 0 0 1/8]},6};
test_fun = @(x) x ./ 6.41;               element11{length(element11)+1}    = {test_fun,{B}, test_fun(B_val)             ,der_from_matrix(1/6.41*(ones(size(B_val)))),n_in_B};
test_fun = @(x) [2 7 6] / x;             element11{length(element11)+1}    = {test_fun,{x}, test_fun(x_val)             ,{[-2/(x_val^2) 0 0 0 0],[-7/(x_val^2) 0 0 0 0],[-6/(x_val^2) 0 0 0 0],[-7/(x_val^2) 0 0 0 0],[-6/(x_val^2) 0 0 0 0]},5};
test_fun = @(x) [2;7;6] / x;             element11{length(element11)+1}    = {test_fun,{x}, test_fun(x_val)             ,{[-2/(x_val^2) 0 0 0 0];[-7/(x_val^2) 0 0 0 0];[-6/(x_val^2) 0 0 0 0]},5};
test_fun = @(x) [2 7 6;8 6 3] / x;       element11{length(element11)+1}    = {test_fun,{x}, test_fun(x_val)             ,{[-2/(x_val^2) 0 0 0 0],[-7/(x_val^2) 0 0 0 0],[-6/(x_val^2) 0 0 0 0];[-8/(x_val^2) 0 0 0 0],[-6/(x_val^2) 0 0 0 0],[-3/(x_val^2) 0 0 0 0]},5};
test_fun = @(x,y) x / y;                 element11{length(element11)+1}    = {test_fun,{x,y},   test_fun(x_val,y_val)   ,{[1/y_val -x_val/(y_val^2) 0 0 0]},5};

%% MATRIX ARITHMETIC OPERATION - TWELTH TEST SERIES MPOWER(^)
test_fun = @(x) x ^ 6;                   element12{1}                      = {test_fun,{x}, test_fun(x_val)             ,{[6.*x_val.^5 0 0 0 0]},5};
test_fun = @(x) x ^ 3;                   element12{length(element12)+1}    = {test_fun,{A}, test_fun(A_val)             ,{[3*a1_val^2+2*a2_val*a3_val, 2*a1_val*a3_val + a3_val*a4_val,2*a1_val*a2_val + a2_val*a4_val, a2_val*a3_val],[2*a3_val*a1_val + a3_val*a4_val, a3_val^2, a4_val^2 + 2*a2_val*a3_val + a1_val^2 + a1_val*a4_val, 2*a3_val*a4_val + a1_val*a3_val];[2*a1_val*a2_val + a2_val*a4_val, a1_val^2 + 2*a2_val*a3_val + a1_val*a4_val + a4_val^2, a2_val^2, a1_val*a2_val + 2*a2_val*a4_val],[a2_val*a3_val, a1_val*a3_val + 2*a3_val*a4_val, a1_val*a2_val + 2*a2_val*a4_val, 3*a4_val^2 + 2*a2_val*a3_val]},4};
test_fun = @(x) 6 ^ x;                   element12{length(element12)+1}    = {test_fun,{x}, test_fun(x_val)             ,{[log(6).*6.^x_val 0 0 0 0]},5};
test_fun = @(x,y) x ^ y;                 element12{length(element12)+1}    = {test_fun,{x,y},   test_fun(x_val,y_val)   ,{[y_val.*x_val.^(y_val-1) log(x_val).*x_val.^y_val 0 0 0]},5};

%% MATRIX ARITHMETIC OPERATION - THIRTEENTH TEST SERIES MTIMES(*)
test_fun = @(x) x * 2;                  element13{1}                       = {test_fun,{x}, test_fun(x_val)          ,{[2 0 0 0 0]},5};
test_fun = @(x) 2 * x;                  element13{length(element13)+1}     = {test_fun,{x}, test_fun(x_val)          ,{[2 0 0 0 0]},5};
test_fun = @(x) x * [2 7 6];            element13{length(element13)+1}     = {test_fun,{x}, test_fun(x_val)          ,{[2 0 0 0 0],[7 0 0 0 0],[6 0 0 0 0]},5};
test_fun = @(x) x * [2;7;6];            element13{length(element13)+1}     = {test_fun,{x}, test_fun(x_val)          ,{[2 0 0 0 0];[7 0 0 0 0];[6 0 0 0 0]},5};
test_fun = @(x) x * [2 7 6;8 6 3];      element13{length(element13)+1}     = {test_fun,{x}, test_fun(x_val)          ,{[2 0 0 0 0],[7 0 0 0 0],[6 0 0 0 0];[8 0 0 0 0],[6 0 0 0 0],[3 0 0 0 0]},5};

test_fun = @(x) x * 8;                  element13{length(element13)+1}     = {test_fun,{b}, test_fun(b_val)          ,{[8 0 0 0 0 0],[0 8 0 0 0 0],[0 0 8 0 0 0],[0 0 0 8 0 0],[0 0 0 0 8 0],[0 0 0 0 0 8]},6};
test_fun = @(x) x * 8;                  element13{length(element13)+1}     = {test_fun,{b'}, test_fun(b_val')        ,{[8 0 0 0 0 0];[0 8 0 0 0 0];[0 0 8 0 0 0];[0 0 0 8 0 0];[0 0 0 0 8 0];[0 0 0 0 0 8]},6};
test_fun = @(x) x * 6.41;               element13{length(element13)+1}     = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix(6.41*(ones(size(B_val)))),n_in_B};

test_fun = @(x) [2 7 6] * x;            element13{length(element13)+1}     = {test_fun,{x}, test_fun(x_val)          ,{[2 0 0 0 0],[7 0 0 0 0],[6 0 0 0 0]},5};
test_fun = @(x) [2;7;6] * x;            element13{length(element13)+1}     = {test_fun,{x}, test_fun(x_val)          ,{[2 0 0 0 0];[7 0 0 0 0];[6 0 0 0 0]},5};
test_fun = @(x) [2 7 6;8 6 3] * x;      element13{length(element13)+1}     = {test_fun,{x}, test_fun(x_val)          ,{[2 0 0 0 0],[7 0 0 0 0],[6 0 0 0 0];[8 0 0 0 0],[6 0 0 0 0],[3 0 0 0 0]},5};

test_fun = @(x) 8 * x;                  element13{length(element13)+1}     = {test_fun,{b}, test_fun(b_val)          ,{[8 0 0 0 0 0],[0 8 0 0 0 0],[0 0 8 0 0 0],[0 0 0 8 0 0],[0 0 0 0 8 0],[0 0 0 0 0 8]},6};
test_fun = @(x) 8 * x;                  element13{length(element13)+1}     = {test_fun,{b'}, test_fun(b_val')        ,{[8 0 0 0 0 0];[0 8 0 0 0 0];[0 0 8 0 0 0];[0 0 0 8 0 0];[0 0 0 0 8 0];[0 0 0 0 0 8]},6};
test_fun = @(x) 6.41 * x;               element13{length(element13)+1}     = {test_fun,{B}, test_fun(B_val)          ,der_from_matrix(6.41*(ones(size(B_val)))),n_in_B};

test_fun = @(x,y) x * y;                element13{length(element13)+1}     = {test_fun,{x,y},   test_fun(x_val,y_val)   ,{[y_val x_val 0 0 0]},5};
test_fun = @(x,y) x * y;                element13{length(element13)+1}     = {test_fun,{x,yy},  test_fun(x_val,yy_val)  ,{[yy_val(1) x_val 0 0 0],[yy_val(2) 0 x_val 0 0]},5};
test_fun = @(x,y) x * y;                element13{length(element13)+1}     = {test_fun,{x,yy'}, test_fun(x_val,yy_val') ,{[yy_val(1) x_val 0 0 0];[yy_val(2) 0 x_val 0 0]},5};
test_fun = @(x,y) x * y;                element13{length(element13)+1}     = {test_fun,{x,B},   test_fun(x_val,B_val)   ,{[x_val+B_val(1,1) 0 0 0 0 0],[B_val(1,2) 0 x_val 0 0 0],[B_val(1,3) 0 0 0 x_val 0];[B_val(2,1) x_val 0 0 0 0],[B_val(2,2) 0 0 x_val 0 0],[B_val(2,3) 0 0 0 0 x_val]},6};

test_fun = @(x,y) x * y;                element13{length(element13)+1}     = {test_fun,{yy,x},  test_fun(yy_val,x_val)  ,{[yy_val(1) x_val 0 0 0],[yy_val(2) 0 x_val 0 0]},5};
test_fun = @(x,y) x * y;                element13{length(element13)+1}     = {test_fun,{yy',x}, test_fun(yy_val',x_val) ,{[yy_val(1) x_val 0 0 0];[yy_val(2) 0 x_val 0 0]},5};
test_fun = @(x,y) x * y;                element13{length(element13)+1}     = {test_fun,{B,x},   test_fun(B_val,x_val)   ,{[x_val+B_val(1,1) 0 0 0 0 0],[B_val(1,2) 0 x_val 0 0 0],[B_val(1,3) 0 0 0 x_val 0];[B_val(2,1) x_val 0 0 0 0],[B_val(2,2) 0 0 x_val 0 0],[B_val(2,3) 0 0 0 0 x_val]},6};

test_fun = @(x) x * [2;7;6;3];          element13{length(element13)+1}     = {test_fun,{c}, test_fun(c_val)           ,{[2 7 6 3]},n_in_c};
test_fun = @(x) x * [2 7 6 3];          element13{length(element13)+1}     = {test_fun,{c'}, test_fun(c_val')         ,{[2 0 0 0],[7 0 0 0],[6 0 0 0],[3 0 0 0];[0 2 0 0],[0 7 0 0],[0 6 0 0],[0 3 0 0];[0 0 2 0],[0 0 7 0],[0 0 6 0],[0 0 3 0];[0 0 0 2],[0 0 0 7],[0 0 0 6],[0 0 0 3]},n_in_c};
test_fun = @(x) x * [2;4];              element13{length(element13)+1}     = {test_fun,{A}, test_fun(A_val)           ,{[2 0 4 0];[0 2 0 4]},4};
test_fun = @(x) x * [2 3;4 1];          element13{length(element13)+1}     = {test_fun,{A}, test_fun(A_val)           ,{[2 0 4 0],[3 0 1 0];[0 2 0 4],[0 3 0 1]},4};

test_fun = @(x) [2;7;6;3] * x;          element13{length(element13)+1}     = {test_fun,{c}, test_fun(c_val)           ,{[2 0 0 0],[0 2 0 0],[0 0 2 0],[0 0 0 2];[7 0 0 0],[0 7 0 0],[0 0 7 0],[0 0 0 7];[6 0 0 0],[0 6 0 0],[0 0 6 0],[0 0 0 6];[3 0 0 0],[0 3 0 0],[0 0 3 0],[0 0 0 3]},n_in_c};
test_fun = @(x) [2 7 6 3] * x;          element13{length(element13)+1}     = {test_fun,{c'}, test_fun(c_val')         ,{[2 7 6 3]},n_in_c};
test_fun = @(x) [2 4] * x;              element13{length(element13)+1}     = {test_fun,{A}, test_fun(A_val)           ,{[2 4 0 0],[0 0 2 4]},4};

test_fun = @(x) [2 3 2 4;4 1 0 3;3 4 6 7; 1 6 4 3]*x;element13{length(element13)+1}={test_fun,{c'}, test_fun(c_val')  ,{[2 3 2 4],[2 3 2 4],[2 3 2 4],[2 3 2 4];[4 1 0 3],[4 1 0 3],[4 1 0 3],[4 1 0 3];[3 4 6 7],[3 4 6 7],[3 4 6 7],[3 4 6 7];[1 6 4 3],[1 6 4 3],[1 6 4 3],[1 6 4 3]},4};
test_fun = @(x) [2 3;4 1] * x;          element13{length(element13)+1}     = {test_fun,{A}, test_fun(A_val)           ,{[2 3 0 0],[0 0 2 3];[4 1 0 0],[0 0 4 1]},4};

%NEW
test_fun = @(x,y) x * y;                element13{length(element13)+1}     = {test_fun,{c,d'}, test_fun(c_val,d_val') ,{[3 1 4 9 4 5 6 7]},8};
test_fun = @(x,y) x * y;                element13{length(element13)+1}     = {test_fun,{c',d}, test_fun(c_val',d_val) ,{[3 0 0 0 4 0 0 0],[1 0 0 0 0 4 0 0],[4 0 0 0 0 0 4 0],[9 0 0 0 0 0 0 4];[0 3 0 0 5 0 0 0],[0 1 0 0 0 5 0 0],[0 4 0 0 0 0 5 0],[0 9 0 0 0 0 0 5];[0 0 3 0 6 0 0 0],[0 0 1 0 0 6 0 0],[0 0 4 0 0 0 6 0],[0 0 9 0 0 0 0 6];[0 0 0 3 7 0 0 0],[0 0 0 1 0 7 0 0],[0 0 0 4 0 0 7 0],[0 0 0 9 0 0 0 7]},8};
test_fun = @(x,y) x * y;                element13{length(element13)+1}     = {test_fun,{e,A},  test_fun(e_val,A_val)  ,{[6 2 0 0 3 4],[0 0 6 2 2 5]},6};
test_fun = @(x,y) x * y;                element13{length(element13)+1}     = {test_fun,{A,e'},  test_fun(A_val,e_val'),{[6 0 2 0 3 2];[0 6 0 2 4 5]},6};

%% MATRIX ARITHMETIC OPERATION - FOURTEENTH TEST SERIES Probability
test_fun = @(x) normcdf(x,0,1);                  element14{1}                       = {test_fun,{x}, test_fun(x_val)  ,{[normpdf(x_val,0,1) 0 0 0 0]},5};
test_fun = @(x) normcdf(x,0,2);                  element14{length(element14)+1}     = {test_fun,{x}, test_fun(x_val)  ,{[normpdf(x_val,0,2) 0 0 0 0]},5};
test_fun = @(x) normcdf(x,-1,2);                 element14{length(element14)+1}     = {test_fun,{x}, test_fun(x_val)  ,{[normpdf(x_val,-1,2) 0 0 0 0]},5};
test_fun = @(x) normcdf(x,pi,2);                 element14{length(element14)+1}     = {test_fun,{x}, test_fun(x_val)  ,{[normpdf(x_val,pi,2) 0 0 0 0]},5};

test_fun = @(x) normcdf(x,0,1);                  element14{length(element14)+1}     = {test_fun,{b}, test_fun(b_val)  ,{[normpdf(b_val(1),0,1)  0 0 0 0 0],[0 normpdf(b_val(2),0,1)  0 0 0 0],[0 0 normpdf(b_val(3),0,1)  0 0 0],[0 0 0 normpdf(b_val(4),0,1)  0 0],[0 0 0 0 normpdf(b_val(5),0,1)  0],[0 0 0 0 0 normpdf(b_val(6),0,1)]},6};
test_fun = @(x) normcdf(x,0,2);                  element14{length(element14)+1}     = {test_fun,{b}, test_fun(b_val)  ,{[normpdf(b_val(1),0,2)  0 0 0 0 0],[0 normpdf(b_val(2),0,2)  0 0 0 0],[0 0 normpdf(b_val(3),0,2)  0 0 0],[0 0 0 normpdf(b_val(4),0,2)  0 0],[0 0 0 0 normpdf(b_val(5),0,2)  0],[0 0 0 0 0 normpdf(b_val(6),0,2)]},6};
test_fun = @(x) normcdf(x,-1,2);                 element14{length(element14)+1}     = {test_fun,{b}, test_fun(b_val)  ,{[normpdf(b_val(1),-1,2) 0 0 0 0 0],[0 normpdf(b_val(2),-1,2) 0 0 0 0],[0 0 normpdf(b_val(3),-1,2) 0 0 0],[0 0 0 normpdf(b_val(4),-1,2) 0 0],[0 0 0 0 normpdf(b_val(5),-1,2) 0],[0 0 0 0 0 normpdf(b_val(6),-1,2)]},6};
test_fun = @(x) normcdf(x,pi,2);                 element14{length(element14)+1}     = {test_fun,{b}, test_fun(b_val)  ,{[normpdf(b_val(1),pi,2) 0 0 0 0 0],[0 normpdf(b_val(2),pi,2) 0 0 0 0],[0 0 normpdf(b_val(3),pi,2) 0 0 0],[0 0 0 normpdf(b_val(4),pi,2) 0 0],[0 0 0 0 normpdf(b_val(5),pi,2) 0],[0 0 0 0 0 normpdf(b_val(6),pi,2)]},6};

test_fun = @(x) normcdf(x,0,1);                  element14{length(element14)+1}     = {test_fun,{b'}, test_fun(b_val'),{[normpdf(b_val(1),0,1)  0 0 0 0 0];[0 normpdf(b_val(2),0,1)  0 0 0 0];[0 0 normpdf(b_val(3),0,1)  0 0 0];[0 0 0 normpdf(b_val(4),0,1)  0 0];[0 0 0 0 normpdf(b_val(5),0,1)  0];[0 0 0 0 0 normpdf(b_val(6),0,1)]},6};
test_fun = @(x) normcdf(x,0,2);                  element14{length(element14)+1}     = {test_fun,{b'}, test_fun(b_val'),{[normpdf(b_val(1),0,2)  0 0 0 0 0];[0 normpdf(b_val(2),0,2)  0 0 0 0];[0 0 normpdf(b_val(3),0,2)  0 0 0];[0 0 0 normpdf(b_val(4),0,2)  0 0];[0 0 0 0 normpdf(b_val(5),0,2)  0];[0 0 0 0 0 normpdf(b_val(6),0,2)]},6};
test_fun = @(x) normcdf(x,-1,2);                 element14{length(element14)+1}     = {test_fun,{b'}, test_fun(b_val'),{[normpdf(b_val(1),-1,2) 0 0 0 0 0];[0 normpdf(b_val(2),-1,2) 0 0 0 0];[0 0 normpdf(b_val(3),-1,2) 0 0 0];[0 0 0 normpdf(b_val(4),-1,2) 0 0];[0 0 0 0 normpdf(b_val(5),-1,2) 0];[0 0 0 0 0 normpdf(b_val(6),-1,2)]},6};
test_fun = @(x) normcdf(x,pi,2);                 element14{length(element14)+1}     = {test_fun,{b'}, test_fun(b_val'),{[normpdf(b_val(1),pi,2) 0 0 0 0 0];[0 normpdf(b_val(2),pi,2) 0 0 0 0];[0 0 normpdf(b_val(3),pi,2) 0 0 0];[0 0 0 normpdf(b_val(4),pi,2) 0 0];[0 0 0 0 normpdf(b_val(5),pi,2) 0];[0 0 0 0 0 normpdf(b_val(6),pi,2)]},6};

%% MATRIX SHAPE      OPERATION - FIRST TEST SERIE
first = 8;          second = 9;         op_element{1}                      = {'two',first,second,I_val(first,second),I_id(first,second)};
first = 12;         second = 14;        op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
first = 2;          second = 4;         op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
first = ':';        second = 4;         op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
first = 6;          second = ':';       op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};

first = ':';        second = ':';       op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
first = 2;          second = 2:5;       op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
first = 4;          second = 2:2:15;    op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
first = 3:7;        second = 6;         op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
first = 3:2:9;      second = 9;         op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};

first = ':';        second = 2:8;       op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
first = ':';        second = 3:4:18;    op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
first = 3:9;        second = ':';       op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
first = 4:6:20;     second = ':';       op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
first = 6:15;       second = 8:12;      op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};

first = 4:6:20;     second = 3:4:18;    op_element{length(op_element)+1}   = {'two',first,second,I_val(first,second),I_id(first,second)};
first = 4;          second = [];        op_element{length(op_element)+1}   = {'one',first,second,I_val(first),I_id(first)};
first = 55;         second = [];        op_element{length(op_element)+1}   = {'one',first,second,I_val(first),I_id(first)};
first = 45:55;      second = [];        op_element{length(op_element)+1}   = {'one',first,second,I_val(first),I_id(first)};
first = 3:5:68;     second = [];        op_element{length(op_element)+1}   = {'one',first,second,I_val(first),I_id(first)};

first = ':';        second = [];        op_element{length(op_element)+1}   = {'one',first,second,I_val(first),I_id(first)};

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
name = {name_1,name_2,name_3,name_4,name_5,name_6,name_7,name_8,name_9,name_10,name_11,name_12,name_13,name_14};
element = {element1,element2,element3,element4,element5,element6,element7,element8,element9,element10,element11,element12,element13,element14};

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
            temp_n_in           = temp_str{5}; 
            reset_tape(init_size_tape);
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
            
            temp_id = temp_eval.id;
            
            for i_row=1:n_row
                for i_col=1:n_col
                    B = reverse_tape(temp_n_in,temp_id(i_row,i_col));
                    %Value
                    if(abs(temp_corr_val(i_row,i_col)-temp_eval.val(i_row,i_col))<10^(-14))
                        val_bool = val_bool*true;
                    else
                        val_bool = val_bool*false;
                    end
                    temp_corr_der = temp_corr_der_str{i_row,i_col};
                    %Derivative
                    if(max(abs(B - temp_corr_der))<10^(-14))
                        der_bool = der_bool*true;
                     else
                        der_bool = der_bool*false;
                    end
                end
            end
            %% OUTPUT
            %Dela upp outputen i delar av fem
            if(mod(i_test-1,5)==0 && i_test~=1)
                fprintf(fid1,'\n');
            end
            if(val_bool)
                if(i_test>=10)
                fprintf(fid1,'The value for test %0.f is \t    correct ',i_test);
                else
                    fprintf(fid1,'The value for test %0.f  is \t    correct ',i_test);
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
                    fprintf(fid1,'and the derivatives for test %0.f is \tNOT        correct.\n',i_test);
                else
                    fprintf(fid1,'and the derivatives for test %0.f  is \tNOT        correct.\n',i_test);
                end
            end
        else
            fprintf(fid1,'*******************%s*********************************************\n',temp_str);
            i_test = 0;
        end
    catch me
        fprintf(fid1,    '*******************ERROR FOR TEST %0.f*************************************************',i_test);
        if(strcmp(me.identifier,'MATLAB:sub2ind:IndexOutOfRange'))
            fprintf(fid1,'CHECK The NUMBER OF INPUT VARIABLES. To solve it extend the init_size_tape in the setup above.'); %This is a known error in the unit test but not in the code
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

for i_test=1:length(op_element)
    
    t_element = op_element{i_test};
    op_case   = t_element{1};
    first     = t_element{2};
    second    = t_element{3};
    r_val     = t_element{4};
    r_id      = t_element{5};
    
    switch op_case
        case 'one'
            b = I(first);
            if(b.val==r_val)
                bool_val = true;
            else
                bool_val = false;
            end
            if(b.id==r_id)
                bool_id = true;
            else
                bool_id = false;
            end
        case 'two'
            b = I(first,second);
            if(b.val==r_val)
                bool_val = true;
            else
                bool_val = false;
            end
            if(b.id==r_id)
                bool_id = true;
            else
                bool_id = false;
            end
    end
    
    if(mod(i_test-1,5)==0 && i_test~=1)
        fprintf(fid2,'\n');
    end
    
    if(bool_val)
        if(i_test>=10)
            fprintf(fid2,'The value for test %0.f is \t    correct ',i_test);
        else
            fprintf(fid2,'The value for test %0.f  is \t    correct ',i_test);
        end
    else
        fprintf(fid2,'The value for test %0.f is \t NOT correct ',i_test);
    end
    
    if(bool_id)
        if(i_test>=10)
            fprintf(fid2,'and the id for test %0.f is \t           correct.\n',i_test);
        else
            fprintf(fid2,'and the id for test %0.f  is \t           correct.\n',i_test);
        end
    else
        if(i_test>=10)
            fprintf(fid2,'and the id for test %0.f is \tNOT        correct.\n',i_test);
        else
            fprintf(fid2,'and the id for test %0.f  is \tNOT        correct.\n',i_test);
        end
    end
    
end

fclose(fid2);



