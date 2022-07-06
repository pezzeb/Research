function [Ko] = f_Ko(Oc,Op)

if(not(isempty(Oc)))
    Oc_str = Oc(:,1);
else
    Oc_str = [];
end

if(not(isempty(Op)))
    Op_str = Op(:,1);
else
    Op_str = [];
end

Ko  = unique([Oc_str;Op_str]);

end