function []   = extendtape_interior(T)
%EXTENDTAPE extends a tape given the extension T
% The difference from extendtape is that the T has a different form here,
% that allows that a row can depend on previous rows in T, which is not the
% case in extendtape. 

global tape
global tape_cur_id

global tape_switch

global tape_Qt
global tape_Qt_cur_id

n_row = size(T,1);
n_col = size(T,2);

if(strcmp(tape_switch,'main_tape'))
    tape = sparse([sparse([tape,sparse(tape_cur_id,n_col-size(tape,2))]);T]);
    tape_cur_id = tape_cur_id + n_row;
elseif(strcmp(tape_switch,'Qt_tape'))
    tape_Qt = sparse([sparse([tape_Qt,sparse(tape_Qt_cur_id,n_col-size(tape_Qt,2))]);T]);
    tape_cur_id = tape_cur_id + n_row;
    tape_Qt_cur_id = tape_Qt_cur_id + n_row;
end

end