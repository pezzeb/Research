function [vector_out] = symmetric_diagonal_mtimes(v_right,id,v1,v2,v3,v4,v5,v6,v7)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if(2>= nargin || nargin>9)
    %ERROR
elseif(1 == id(1))
    vector_out = v1.*v_right;
else
    vector_out = zeros(size(v_right));
end
id_length = length(id);

if(id_length>=2)
    vector_out = f(v_right,id(2),v2,vector_out);
end
if(id_length>=3)
    vector_out = f(v_right,id(3),v3,vector_out);
end
if(id_length>=4)
    vector_out = f(v_right,id(4),v4,vector_out);
end
if(id_length>=5)
    vector_out = f(v_right,id(5),v5,vector_out);
end
if(id_length>=6)
    vector_out = f(v_right,id(6),v6,vector_out);
end
if(id_length>=7)
    vector_out = f(v_right,id(7),v7,vector_out);
end
end



function [out] = f(r,id,v,res)
res(id:end) = res(id:end) + r(1:end-id+1).*v;
res(1:end-id+1) = res(1:end-id+1) + r(id:end).*v;
out = res;
end