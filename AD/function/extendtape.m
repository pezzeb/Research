function []   = extendtape(T,ext_type)
%EXTENDTAPE extends a tape given the extension T
% All rows can only depend on intermediate variable that are defined in the
% existing tape and not in the extension T. 

if(1==nargin)
    ext_type = 'normal';
end

global tape
global tape_cur_id


n_row = size(T,1);
n_col = size(T,2);
if(strcmp(ext_type,'normal'))
    if(tape_cur_id==n_col)
        tape_cur_id = tape_cur_id + n_row;
        tape = sparse([sparse([tape; T]),sparse(tape_cur_id,n_row)]);
    else
        error('ERROR')
    end
elseif(strcmp(ext_type,'ext_Qt'))
    tape = sparse([tape,sparse(tape_cur_id,n_row);T,sparse(n_row,tape_cur_id+n_row-n_col)]);
    tape_cur_id = tape_cur_id + n_row;
elseif(strcmp(ext_type,'interpol_pricing'))
    tape = sparse([[tape;T],sparse(tape_cur_id,size(T,1))]);
end
end