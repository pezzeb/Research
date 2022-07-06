function [] = reset_tape(n_in)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global tape;
global tape_cur_id;

global tape_switch

global tape_Qt
global tape_Qt_cur_id


if(strcmp(tape_switch,'main_tape'))
    if(1==nargin)
        tape = sparse(n_in,n_in);
        tape_cur_id = n_in;
    else
        tape            = sparse(n_in,n_in);
        tape_cur_id     = n_in;
    end
elseif(strcmp(tape_switch,'Qt_tape'))
    tape_Qt = sparse(0,0);
    tape_Qt_cur_id = 0;
end
end

