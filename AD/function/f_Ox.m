function [Oc_out,Op_out,S_ini_out] = f_Ox(type_in,excel_path,Oc_in,Op_in,S_ini_in)

if(strcmp('fix',type_in))
    Oc_out    = Oc_in;
    Op_out    = Op_in;
    
    S_ini_out = S_ini_in;
else(strcmp('load',type_in));
    load(excel_path);
    Oc_out    = Oc;
    Op_out    = Op;
    
    Oc_out(:,2)    = Oc_out(:,2)/252;
    Op_out(:,2)    = Op_out(:,2)/252;
    
    S_ini_out = S_ini;
end


end