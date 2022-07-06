function [output_var] = reverse_tape(n_in_var,id)
global tape;


if(id>size(tape,2))
        tape_fn = tape(1:id,:);
else
    tape_fn = tape(1:id,1:id);
end


tape_fn = [tape_fn,sparse(id,id-size(tape_fn,2))];
tape_fn = tape_fn - speye(id,id);

n_rhs = size(tape_fn,1);    
rhs = zeros(1,n_rhs);
rhs(end) = -1;

reversed_tape = rhs/tape_fn;

output_var = reversed_tape(1:n_in_var);
end