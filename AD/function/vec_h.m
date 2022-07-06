function [out] = vec_h(M)
%VEC_V Horizontal vectorization
%   Does a horizontal vectorization of the matrix M.
%   Returns a column vector.

A = M';
out = A(:);

end