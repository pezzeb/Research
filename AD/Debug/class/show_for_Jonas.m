x = forADm(2,[1,0,0]);
y = forADm(3,[0,1,0]);
z = forADm(5,[0,0,1]);
% 
x_mult_y = x*y
x_plus_y = x+y
x_minus_y = x-y
x_mult_x = x*x
x_pow_3 = x^3
x_pow_z = x^z

h = sin(x_pow_z)
x_exp = exp(x)
x_sin = sin(x)
x_cos = cos(x)

x_ny = x_mult_y %Visa att det går att tilldela

%
xx = forADm([5;6],[1,0,0;0,1,0]);
yy = forADm([3;2],[0,1,0;0,1,0]);
zz = forADm(5,[0,0,1]);

omega = forADm([3;4;5],eye(3))
%kod 
%loop
omega = forADm(omega.val,eye(3));

mats = [1 1 1;2 2 2;3 3 3]*omega

%%
vector_x = [x;y] %Så här är det INTE möjligt att göra med min kod! Ett objekt skapas men det det går inte att använda sedan.

