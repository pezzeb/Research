classdef revADm
    %revADm is the class for autodifferentiation
    %The class is used as objects when the reverse autodifferentiations is
    %used
    properties
        val         %The value of the object
        id          %The id of the object
        %All instances have acsess to a global tape.
    end
    methods
        
        %% CONSTRUCTOR & VALUE
        function obj = revADm(a,b)
            
            if nargin ==1
                obj.val = a;
                obj.id  = zeros(size(a));
            elseif  nargin == 2
                
                if(size(a) == size(b))
                    obj.val = a;
                    obj.id = b;
                else
                    error('USER-ERROR: Not correct dimension in revADm constructor');
                end
            else
                error('USER-ERROR: Not correct number of arguments in AD constructor');
            end
        end
        function vec = double(obj)
            %Convert revAD object to matrix of doubles.
            vec = [ obj.val, obj.id ];
        end
        
        %% One dimensional operators
        %These are functions that only take one argument.
        function h = uminus(u)
            new_id = one_arg_ins(u.id,-ones(size(u.val)));
            h = revADm(-(u.val),new_id);
        end
        function h = exp(u)   
            new_id = one_arg_ins(u.id,exp(u.val));
            h = revADm(exp(u.val),new_id);
        end
        function h = log(u)   
            new_id = one_arg_ins(u.id,1./(u.val));
            h = revADm(log(u.val),new_id);
        end
        function h = sin(u)   
            new_id = one_arg_ins(u.id,cos(u.val));
            h = revADm(sin(u.val),new_id);
        end
        function h = sinh(u)  
            new_id = one_arg_ins(u.id,cosh(u.val));
            h = revADm(sinh(u.val),new_id);
        end
        function h = asin(u)  
            new_id = one_arg_ins(u.id,1./sqrt(1-u.val.^2));
            h = revADm(asin(u.val),new_id);
        end
        function h = asinh(u) 
            new_id = one_arg_ins(u.id,1./sqrt(u.val.^2+1));
            h = revADm(asinh(u.val),new_id);
        end
        function h = cos(u)   
            new_id = one_arg_ins(u.id,-sin(u.val));
            h = revADm(cos(u.val),new_id);
        end
        function h = cosh(u)  
            new_id = one_arg_ins(u.id,sinh(u.val));
            h = revADm(cosh(u.val),new_id);
        end
        function h = acos(u)  
            new_id = one_arg_ins(u.id,-1./sqrt(1-u.val.^2));
            h = revADm(acos(u.val),new_id);
        end
        function h = acosh(u) 
            new_id = one_arg_ins(u.id,1./sqrt(u.val.^2-1));
            h = revADm(acosh(u.val),new_id);
        end
        function h = tan(u)   
            new_id = one_arg_ins(u.id,sec(u.val).^2);
            h = revADm(tan(u.val),new_id);
        end
        function h = tanh(u)  
            new_id = one_arg_ins(u.id,1 - tanh(u.val).^2);
            h = revADm(tanh(u.val),new_id);
        end
        function h = atan(u)  
            new_id = one_arg_ins(u.id,1./(1+u.val.^2));
            h = revADm(atan(u.val),new_id);
        end
        function h = atanh(u) 
            new_id = one_arg_ins(u.id,1./(1-u.val.^2));
            h = revADm(atanh(u.val),new_id);
        end
        function h = csc(u)   
            new_id = one_arg_ins(u.id,-cot(u.val).*csc(u.val));
            h = revADm(csc(u.val),new_id);
        end
        function h = csch(u)  
            new_id = one_arg_ins(u.id,-coth(u.val).*csch(u.val));
            h = revADm(csch(u.val),new_id);
        end
        function h = acsc(u)  
            new_id = one_arg_ins(u.id,-1./(abs(u.val).*sqrt(u.val.^2-1)));
            h = revADm(acsc(u.val),new_id);
        end
        function h = acsch(u) 
            new_id = one_arg_ins(u.id,-1./(abs(u.val).*sqrt(1+u.val.^2)) );
            h = revADm(acsch(u.val),new_id);
        end
        function h = cot(u)   
            new_id = one_arg_ins(u.id,-csc(u.val).^2);
            h = revADm(cot(u.val),new_id);
        end
        function h = coth(u)  
            new_id = one_arg_ins(u.id,1-coth(u.val).^2);
            h = revADm(coth(u.val),new_id);
        end
        function h = acot(u)  
            new_id = one_arg_ins(u.id,-1./(1+u.val.^2));
            h = revADm(acot(u.val),new_id);
        end
        function h = acoth(u) 
            new_id = one_arg_ins(u.id,1./(1-u.val.^2));
            h = revADm(acoth(u.val),new_id);
        end
        function h = sec(u)   
            new_id = one_arg_ins(u.id,sec(u.val).*tan(u.val));
            h = revADm(sec(u.val),new_id);
        end
        function h = sech(u)  
            new_id = one_arg_ins(u.id,-tanh(u.val).*sech(u.val));
            h = revADm(sech(u.val),new_id);
        end
        function h = asec(u)  
            new_id = one_arg_ins(u.id,1./(abs(u.val).*sqrt(u.val.^2-1)));
            h = revADm(asec(u.val),new_id);
        end
        function h = asech(u) 
            new_id = one_arg_ins(u.id,-1./(abs(u.val).*sqrt(1-u.val.^2)));
            h = revADm(asech(u.val),new_id);
        end
        function h = sqrt(u)
            new_id = one_arg_ins(u.id,1./(2.*sqrt(u.val)));
            h = revADm(sqrt(u.val),new_id);
        end
        function h = gamma(u)
            tmp    = (gamma(u.val+0.001)-gamma(u.val-0.001))/0.002;
            tmpnew = gamder(u.val);
            new_id = one_arg_ins(u.id,tmp); %tmp <- gamma(u.val).*psi(u.val)
            h = revADm(gamma(u.val),new_id);
        end
        function h = normcdf(u,mu,sigma) %Lite "fuskat" eftersom att man måste ange mu och sigma
            new_id = one_arg_ins(u.id,normpdf(u.val,mu,sigma));
            h = revADm(normcdf(u.val,mu,sigma),new_id);
        end
        function h = conj(z)
            new_id = intermediatereverse(z);
            h = revADm(conj(z.val),new_id);
        end
        function h = real(z)
            w = conj(z);
            h = (z + w)/2;
        end
        function h = imag(z)
            h = (z - conj(z))/(2*1i);
        end
        
        function new_id = one_arg_temp(u)
            global tape;
            global tape_cur_id;
            T = sparse(numel(u.id),tape_cur_id);
            T = conj(tape(u.id,:));
            new_id = tape_cur_id+1:1:tape_cur_id+numel(u.id);
            extendtape(T);
        end
        function new_id = intermediatereverse(u)
            global tape_cur_id;
            global n_var_revADm;
            uid = u.id;
            numelu = numel(uid);
            T = sparse(numelu,tape_cur_id);
            temprev = sparse(numelu,n_var_revADm);
            
            for i=1:numelu
                temprev(i,:) = reverse_tape(n_var_revADm,uid(i)); 
            end
            
            T(1:numelu,1:n_var_revADm) = conj(temprev);
            new_id = (tape_cur_id+1:1:tape_cur_id+numel(uid))'; %Antar att man får in en kolumnvektor tror jag
            
            extendtape(T);
            
        end
        
        
        %% Two dimensional operators of type interior
        function h = rdivide(u,v)
            [~,~,~,~,u_val] = var_type(u);
            [~,~,~,~,v_val] = var_type(v);
            
            %Calculate the (eventual) derivatives
            T1_val = 1./v_val;
            T2_val = -u_val./(v_val.^2);
            %Calcualtes the value
            h_val  = u_val./v_val;
            
            h=two_arg_interior(u,v,h_val,T1_val,T2_val);
        end
        function h = mrdivide(u,v)
            
            [vtype,~,~,~,~] = var_type(v);
            if(matrixType.Scalar == vtype)
                
                h = u./ v;
            else
                error('USER-ERROR: revADm mrdivide, not implemented functionallity');
            end
        end
        function h = mldivide(u,v)
            
            [n,m] = size(u);
            if(n==m)
                L = own_chol(u);
                u1 = lower_trig_solverSpeed(L,v);
                h  = upper_trig_solverSpeed(L',u1);
            else
                error('USER-ERROR: The matrix in mldivide is not square and it is not implemented');
            end
            
        end
        function h = plus(u,v)
            [~,~,~,~,u_val] = var_type(u);
            [~,~,~,~,v_val] = var_type(v);
            
            %Calculate the (eventual) derivatives
            T1_val = ones(size(u_val));
            T2_val = ones(size(v_val));
            %Calcualtes the value
            h_val  = u_val + v_val;
            
            h=two_arg_interior(u,v,h_val,T1_val,T2_val);
        end
        function h = minus(u,v)
            h = u + (-v);
        end
        function h = times(u,v)
            [~,~,~,~,u_val] = var_type(u);
            [~,~,~,~,v_val] = var_type(v);
            %Calculate the (eventual) derivatives
            T1_val = v_val.*ones(size(u_val));
            T2_val = u_val.*ones(size(v_val));
            %Calcualtes the value
            h_val  = u_val.*v_val;
            
            h=two_arg_interior(u,v,h_val,T1_val,T2_val);
            
        end
        function h = power(u,v)
            [~,~,~,~,u_val] = var_type(u);
            [~,~,~,~,v_val] = var_type(v);
            
            %Calculate the (eventual) derivatives
            T1_val = v_val.*u_val.^(v_val-1);
            T2_val = log(u_val).*u_val.^v_val;
            %Calcualtes the value
            h_val  = u_val.^v_val;
            
            h=two_arg_interior(u,v,h_val,T1_val,T2_val);
        end
        
        %% Two dimensional operators of special type
        function h = mpower(u,v)
            [utype,~,~,~,~] = var_type(u);
            [vtype,is_v,~,~,v_val] = var_type(v);
            
            if(matrixType.Scalar == utype && matrixType.Scalar == vtype)
                %if both is scalar use an already define function
                h = u.^v;
            elseif(matrixType.Scalar == vtype && ~is_v)
                if(mod(v_val,1)==0)
                    %multiply use of the defined function *
                    prod_org = u;
                    prod     = u;
                    for i_pow=1:v-1
                        prod = prod*prod_org;
                    end
                    h = prod;
                else
                    error('USER-ERROR not integer exponent');
                end
            else
                error('Not implemented functionallity of mpower');
            end
        end
        function h = mtimes(u,v)
            %MTIMES is a special function than all the other
            global tape_cur_id;
            
            [utype,is_u,n_row_u,n_col_u,u_val] = var_type(u);
            [vtype,is_v,n_row_v,n_col_v,v_val] = var_type(v);
            
            %Check if any index is zero.
            all_id_bool = true;
            if(is_v)
                if(not(nnz(v.id)==numel(v.id)))
                    all_id_bool = false;
                end
            end
            if(is_u)
                if(not(nnz(u.id)==numel(u.id)))
                    all_id_bool = false;
                end
            end
            
            if(matrixType.Scalar==utype)
                %reuse code
                h = u.*v;
            elseif(matrixType.Scalar==vtype)
                %reuse code
                h = u.*v;
            elseif(n_col_u == n_row_v && all_id_bool)
                if(is_v && is_u)
                    T1 = with_v(n_row_u, n_col_u, n_row_v, n_col_v, u, v, u_val, v_val);
                    T2 = with_u(n_row_u, n_col_u, n_row_v, n_col_v, u, v, u_val, v_val);
                    T = T1 + T2;
                elseif(is_u)
                    T = with_v(n_row_u, n_col_u, n_row_v, n_col_v, u, v, u_val, v_val);
                elseif(is_v)
                    T = with_u(n_row_u, n_col_u, n_row_v, n_col_v, u, v, u_val, v_val);
                else
                    error('ERROR');
                end
                
                h_val = u_val*v_val;
                
                n_new = length(h_val(:));
                new_id_t = tape_cur_id+1:1:tape_cur_id + n_new;
                new_id = reshape(new_id_t,n_row_u,n_col_v);     %reshape to match
                
                h = revADm(h_val,new_id);
                extendtape(T);
            elseif(n_col_u == n_row_v && matrixType.ColumnVector==vtype)
                h = sum(u'.*repmat(v,1,size(u,1)),1)';
                
            elseif(n_col_u == n_row_v && matrixType.RowVector==vtype)
                h = sum(repmat(u',1,size(v,2)).*v,1);
                
            else
                error('ERROR');
            end
        end
        
        function h = speed_mtimes_sparse(u,v)
            global tape_cur_id
            v_val = v.val;
            h_val = u*v_val;
            
            n_id = numel(h_val);
            
            T = sparse(n_id,tape_cur_id);
            T(1:n_id,1:n_id) = u;
            
            ht = revADm(h_val,(tape_cur_id+1:1:tape_cur_id+n_id)');
            h  = ht';
            extendtape(T);
            
        end
        
        %% SPECIAL FUNCTION
        function y = normpdf(x,mu,sigma)
            
            %FROM MATLAB DEFINED FUNCTION
            if nargin<1
                error(message('stats:normpdf:TooFewInputs'));
            end
            if nargin < 2
                mu = 0;
            end
            if nargin < 3
                sigma = 1;
            end
            
            try
                y = exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);
            catch
                error(message('stats:normpdf:InputSizeMismatch'));
            end
        end
        function Y_rep = repmat(X,m,n)
            Y_rep = revADm(repmat(X.val,m,n),repmat(X.id,m,n));
        end
        
        function y = sum(x,dim)
            
            global tape_cur_id;
            
            if nargin==1
                dim = 1;
            end
            
            %FALLET MED NOLLOR ÄR EJ IMPLEMENTERAT
            
            if(1==dim)
                
                n_min     = max(x.id,[],1);
                n_min_temp= n_min;
                n_min_temp(n_min>0)=1;
                n_min_mat = logical(repmat(n_min_temp,size(x.id,1),1));
                
                single_row = 1:1:nnz(n_min);
                row_mat    = repmat(single_row,size(x.id,1),1);
                
                x_id_t = reshape(x.id(n_min_mat),size(row_mat,1),size(row_mat,2));
                
                n_el   = nnz(n_min);
                
                
                new_id = tape_cur_id+1:1:tape_cur_id+n_el;
                
                row_id = row_mat(x_id_t>0);
                col_id = x.id(x.id>0);
                T      = sparse(n_el,tape_cur_id);
                idx    = sub2ind(size(T),row_id,col_id);
                T(idx) = 1;
                
                extendtape(T);
                new_id_ext = zeros(size(n_min_temp));
                new_id_ext(logical(n_min_temp)) = new_id;
                y = revADm(sum(x.val,dim),new_id_ext);
                
            elseif(2==dim)
                y = sum(x',1)';
            else
                error('USER-ERROR not implemented in SUM revADm');
            end
        end
        
        function y = spdiags(M,pl,n_row,n_col)
            y = revADm(spdiags(M.val,pl,n_row,n_col),spdiags(M.id,pl,n_row,n_col));
        end
        function y = diag(M)
            y = revADm(diag(M.val),diag(M.id));
        end
        
        %% TRANSPOSE and SHAPING
        function h = ctranspose(u)
            h = revADm(u.val',u.id');
        end
        function h = transpose(u)
            h = revADm(u.val.',u.id.');
        end
        function [r,c] = size(u,dim)
            if(exist('dim','var'))
                r = size(u.val,dim);
            else
                [r,c] = size(u.val);
            end
        end
        function [l]   = length(v)
            l = length(v.val);
        end
        function [re]  = reshape(x,m,n)
            re = revADm(reshape(x.val,m,n), reshape(x.id,m,n));
        end
        
        %% LOCIGICAL OPERATIONS
        function h = le(u,v)
            [~,is_u,~,~,~] = var_type(u);
            [~,is_v,~,~,~] = var_type(v);
            
            if(is_u && is_v)
                h = u.val>=v.val;
            elseif(is_u)
                h = u.val>=v;
            elseif(is_v)
                h= u>=v.val;
            else
                error('USER-ERROR in revADm le');
            end
        end
        function hlogi = eq(u,v)
            [~,is_u,~,~,~] = var_type(u);
            [~,is_v,~,~,~] = var_type(v);
            
            if(is_u && is_v)
                hlogi = u.val==v.val;
            elseif(is_u)
                hlogi = u.val==v;
            elseif(is_v)
                hlogi = u==v.val;
            else
                error('USER-ERROR in revADm le');
            end
            
        end
        
        %% SPEED FUNCTIONS
        function h = speed_func(x,y,z,w,q)
            global tape_cur_id;
            
            new_id = tape_cur_id+1;
            
            T = [4 6 16 1 1];
            extendtape(T);
            
            h = revADm(4*x.val+6*y.val+16*z.val+w.val+q.val,new_id);
        end
        
        function h = normpdfSPEED(x,mu,sigma)
            global tape_cur_id
            %FROM MATLAB DEFINED FUNCTION
            if nargin<1
                error(message('stats:normpdf:TooFewInputs'));
            end
            if nargin < 2
                mu = 0;
            end
            if nargin < 3
                sigma = 1;
            end
            
            try
                n = length(x);
                T = sparse(n,tape_cur_id);
                new_id = tape_cur_id + (1:1:n)';
                
                T(:,sigma.id) = 1/sigma.val.*(((x-mu.val)./sigma.val).^2-1).*normpdf(x,mu.val,sigma.val);
                T(:,mu.id)    = (x-mu.val)/(sigma.val.^2).*normpdf(x,mu.val,sigma.val);
                extendtape(T);
                h = revADm(normpdf(x,mu.val,sigma.val),new_id);
            catch me
                error(message('stats:normpdf:InputSizeMismatch'));
            end
        end
        function mc = normal_central_moments(sigma,n_m)
            global tape_cur_id;
            switch n_m
                case 1
                    mc = revADm(0,0);
                case 2
                    new_id  = [0;tape_cur_id+1];
                    
                    T           = sparse(1,tape_cur_id);
                    T(sigma.id) = 2*sigma.val;
                    extendtape(T);
                    
                    mc = revADm([0;sigma.val^2],new_id);
                case 3
                    new_id  = [0;tape_cur_id+1;0];
                    
                    T           = sparse(1,tape_cur_id);
                    T(sigma.id) = 2*sigma.val;
                    extendtape(T);
                    
                    mc = revADm([0;sigma.val^2;0],new_id);
                case 4
                    new_id  = [0;tape_cur_id+1;0;tape_cur_id+2];
                    
                    T           = sparse(2,tape_cur_id);
                    T(:,sigma.id) = [2*sigma.val;12*sigma.val^3];
                    extendtape(T);
                    
                    mc = revADm([0;sigma.val^2;0;3.*sigma.val^4],new_id);
                case 5
                    new_id  = [0;tape_cur_id+1;0;tape_cur_id+2;0];
                    
                    T           = sparse(2,tape_cur_id);
                    T(:,sigma.id) = [2*sigma.val;12*sigma.val^3];
                    extendtape(T);
                    
                    mc = revADm([0;sigma.val^2;0;3.*sigma.val^4;0],new_id);
                case 6
                    new_id  = [0;tape_cur_id+1;0;tape_cur_id+2;0;tape_cur_id+3];
                    
                    T           = sparse(3,tape_cur_id);
                    T(:,sigma.id) = [2*sigma.val;12*sigma.val^3;90*sigma.val^5];
                    extendtape(T);
                    
                    mc = revADm([0;sigma.val^2;0;3.*sigma.val^4;0;15.*sigma.val^6],new_id);
                otherwise
                    %ERROR
            end
        end
        
        %% COLON & END
        %https://www.mathworks.com/matlabcentral/newsreader/view_thread/257044
        %http://www.mathworks.com/matlabcentral/newsreader/view_thread/19174
        
        %ÖVERLAGRA INTE NUMEL TÅ FÅR MAN PROBLEM, EFTERSOM ATT DEN ANVÄNDS NÄR MATLAB SJÄLV SKA KOLLA OM MAN HAR FÖR MÅNGA OUTPUT ARGUMENT
        
        
        
        % LÄNKAR TILL END
        %http://compgroups.net/comp.soft-sys.matlab/using-subsref-subsasgn-and-substruct/961250
        %http://www.dsprelated.com/groups/matlab/show/3650.php
        
        %BRA LÄNKAR TILL SUBSREF
        %https://www.google.se/webhp?sourceid=chrome-instant&ion=1&espv=2&ie=UTF-8&safe=active&ssui=on#q=overload+colon+matlab&start=20&safe=active&ssui=on
        %http://se.mathworks.com/help/matlab/ref/colon.html
        %https://www.google.se/webhp?sourceid=chrome-instant&ion=1&espv=2&ie=UTF-8&safe=active&ssui=on#q=different+colons+in+matlab+overload&safe=active&ssui=on
        %http://www.dsprelated.com/groups/matlab/show/3650.php
        %https://www.google.se/webhp?sourceid=chrome-instant&ion=1&espv=2&ie=UTF-8&safe=active&ssui=on#q=substruct+matlab&start=10&safe=active&ssui=on
        %http://compgroups.net/comp.soft-sys.matlab/using-subsref-subsasgn-and-substruct/961250
        %https://www.mathworks.com/matlabcentral/newsreader/view_thread/257044
        %https://www.google.se/webhp?sourceid=chrome-instant&ion=1&espv=2&ie=UTF-8&safe=active&ssui=on#q=subsref+overload+end&safe=active&ssui=on
        %http://stackoverflow.com/questions/20863050/why-does-matlab-throw-a-too-many-output-arguments-error-when-i-overload-subsre
        %http://www.mathworks.com/matlabcentral/newsreader/view_thread/19174
        
        %% SUBSREF & SUBSASGN
        function h = subsref(u,s)
            %SUBSREF is the ., (), .() in essence the colon operator
            %The IN-parameters:
            %u  :   The object
            %s  :   The operator
            %The OUT-parameters
            %h  :   The value to return
            
            switch_type = [s.type];
            switch_subs = [s.subs];
            
            switch switch_type
                case '()'   %obj(??)
                    if(length(s.subs)==1)
                        h = revADm(u.val(switch_subs{1}),u.id(switch_subs{1})); %Colon or linearindex
                    else
                        h = revADm(u.val(switch_subs{1,1},switch_subs{1,2}),u.id(switch_subs{1,1},switch_subs{1,2})); %all other for cases of the type (.,.)
                    end     %obj.??
                case'.'
                    switch s.subs
                        case 'val'
                            h=u.val;
                        case 'id'
                            h=u.id;
                        otherwise
                            error('no such property');
                    end
                case'.()'   %obj.data()
                    switch_subs_val_id = switch_subs{1};
                    switch_subs_rest   = [switch_subs{2:end}];
                    
                    %Determine if val or id (or error)
                    if(strcmp(switch_subs_val_id,'val'))
                        val_id = u.val;
                    elseif(strcmp(switch_subs_val_id,'id'))
                        val_id = u.id;
                    else
                        error('error');
                    end
                    
                    %apply the operation inside the paranthesis
                    if(length(switch_subs_rest)==1)
                        h = val_id(switch_subs_rest);
                    else
                        h = val_id(switch_subs_rest(1,1),switch_subs_rest(1,2));
                    end
                otherwise
                    error('error');
            end
        end
        function h = subsasgn(u,s,b)
            %SUBSASGN is the same as SUBSREF but with a right hand side
            %Can only be used if the right handside writes information to
            %ids in the left hand side that have the current id 0;
            switch_type = [s.type];
            switch_subs = [s.subs];
            
            switch switch_type
                case '()'
                    if(length(s.subs)==1)
                        if(0==(u.id(switch_subs{1})))
                            valt = u.val;
                            valt(switch_subs{1}) = b.val;
                            
                            idt  = u.id;
                            idt(switch_subs{1})  = b.id;
                            
                            h = revADm(valt,idt);
                        else
                            error('USER-DEFINE in subsasgn, already define')
                        end
                    else
                        if(0==(u.id(switch_subs{1},switch_subs{2})))
                            valt = u.val;
                            valt(switch_subs{1},switch_subs{2}) = b.val;
                            
                            idt  = u.id;
                            idt(switch_subs{1},switch_subs{2}) = b.id;
                            
                            h = revADm(valt,idt);
                        else
                            error('USER-DEFINE in subsasgn, already define');
                        end
                    end
                otherwise
                    error('error subsasgn');
            end
        end
    end
end

%HELP FUNCTIONS TO MTIMES
function T = with_v(n_row_u, n_col_u, n_row_v, n_col_v, u, v, u_val, v_val)
% WITH_V is the help function in the operation u*v for the case when v is
% constant
global tape;
global tape_cur_id;

n_el = n_row_u*n_col_v;
T = zeros(n_el,tape_cur_id);

for i_r=1:n_el
    row_id((i_r-1)*n_col_u+1:i_r*n_col_u,1) = i_r;
end

for j_c=1:n_row_u
    col_id_t((j_c-1)*n_col_u+1:j_c*n_col_u,1) = u.id(j_c,:)';
end

col_id = repmat(col_id_t,n_col_v,1);
idx = sub2ind(size(T),row_id,col_id);

for k=1:n_col_v
    temp  = repmat(v_val(:,k),n_row_u,1);
    T_der((k-1)*size(temp,1)+1:k*size(temp,1),1) = temp;
end

T(idx) = T_der;
end
function T = with_u(n_row_u, n_col_u, n_row_v, n_col_v, u, v, u_val, v_val)
% WITH_U is the help function in the operation u*v for the case when u is
% constant

global tape;
global tape_cur_id;

n_el = n_row_u*n_col_v;
T = zeros(n_el,tape_cur_id);

for i_r=1:n_el
    row_id((i_r-1)*n_col_u+1:i_r*n_col_u,1) = i_r;
end

n_t = n_row_u*n_col_u;

for j_c=1:n_col_v
    temp = repmat(v.id(:,j_c),n_row_u,1);
    col_id((j_c-1)*n_t + 1:j_c*n_t,1) = temp;
end

idx = sub2ind(size(T),row_id,col_id);
temp_der = vec_h(u_val);
T_der = repmat(temp_der,n_col_v,1);
T(idx) = T_der;

end

%Help function to MLDIVIDE
function [h] = own_chol(A)
global tape_cur_id
A_val = A.val;

L = chol(A_val,'lower');
sL = size(L,2);

%new_id
temp    = 1:1:(sL^2+sL)/2;
N = tril(ones(sL,sL));
N(N>0) = temp+tape_cur_id;

new_id  = N;

T_L= zeros((sL^2+sL)/2,(sL^2+sL)/2);
start_row   = 1;
stop_row    = sL;
inner_pos   =1;
current_col =1;

for outer_col=sL:-1:1
    v=-L(inner_pos:end,1:inner_pos)./L(inner_pos,inner_pos);
    
    if(not(start_row==stop_row))
        old_off = 0;
        I = eye(outer_col);
        for j=1:inner_pos-1
            
            T_L(start_row+1:stop_row,current_col+old_off) = v(2:end,j);
            
            T_L(start_row:stop_row,current_col+old_off:current_col+old_off+outer_col - 1) =...
                T_L(start_row:stop_row,current_col+old_off:current_col+old_off+outer_col - 1) + I*v(1,j);
            old_off = old_off + (sL-j);
        end
        T_L(start_row+1:stop_row,current_col+old_off) = v(2:end,inner_pos);
    else
        off = 0;
        for jj=1:inner_pos-1
            T_L(end,current_col+off) = v(jj);
            off = off + sL-jj;
        end
    end
    inner_pos = inner_pos+1;
    current_col = current_col+1;
    start_row = start_row + outer_col;
    stop_row  = stop_row  + outer_col-1;
end

t_A= zeros(1,(sL^2+sL)/2);
start_pos = 1;
stop_pos  = sL;
L_pos = 1;
for outer_col=sL:-1:1
    t_A(start_pos:stop_pos) = [1/2 ones(1,outer_col-1)]/L(L_pos,L_pos);
    
    start_pos = start_pos + outer_col;
    stop_pos  = stop_pos  + outer_col-1;
    L_pos = L_pos + 1;
    
end

T_A = diag(t_A);

Tt = sparse(nnz(new_id),tape_cur_id);

At = tril(A.id);

Tt(:,At(At>0)) = T_A;
T  = sparse([Tt,T_L]);
extendtape_interior(T);

h = revADm(L,new_id);

end
function [h] = lower_trig_solverSpeed(A,b)
global tape_cur_id
A_val = A.val;
A_id  = A.id;
a_id  = A_id(A_id>0);

b_val = b.val;
b_id  = b.id;

v = A_val\b_val;

r = diag(A_val);

n_A_row = size(A,1);
len_v   = length(v);

R = -A_val./repmat(r,1,n_A_row);
R = R + 2*eye(n_A_row);

T_b = diag(1./r);

for k=2:length(v);
    for j=1:k-1
        T_b(k,j) = R(k,1:k-1)*T_b(1:k-1,j);
    end
end

nnz_el = n_A_row+(n_A_row^2-n_A_row)/2;

T_A = sparse(len_v,nnz_el);

start = 1;
stop  = n_A_row;

for kk=1:len_v
    T_A(:,start:stop) = T_b(:,kk:end)*(-v(kk));
    
    start = start + (n_A_row+1-kk);
    stop  = stop  + (n_A_row-kk);
end


T = sparse(len_v,tape_cur_id);

T(:,logical(b_id)) = T_b(:,1:end-1);
T(:,logical(a_id)) = T(:,logical(a_id)) + T_A;

new_id = tape_cur_id + (1:1:len_v)';

extendtape(T);
h = revADm(v,new_id);

end
function [h] = upper_trig_solverSpeed(U,b)
global tape_cur_id
U_val = U.val;
U_id  = U.id;
u_id  = U_id(U_id>0);

b_val = b.val;
b_id  = b.id;

v = U_val\b_val;

r = diag(U_val);

n_U_row = size(U,1);
len_v   = length(v);

R = -U_val./repmat(r,1,n_U_row);
R = R + 2*eye(n_U_row);

T_b = diag(1./r);

for j=length(v):-1:2;
    for k=j-1:-1:1
        T_b(k,j) = R(k,k+1:j)*T_b(k:j-1,j);
    end
end

nnz_el = n_U_row+(n_U_row^2-n_U_row)/2;

T_U = sparse(len_v,nnz_el);

start_col = 1;
stop_col  = 1;

for jj=1:len_v
    T_U(:,start_col:stop_col) = T_b(:,1:jj)*(-v(jj));
    
    start_col = start_col + jj;
    stop_col  = stop_col  + jj+1;
end

T = sparse(len_v,tape_cur_id);

T(:,logical(b_id)) = T_b;
T(:,logical(u_id)) = T(:,logical(u_id)) + T_U;

new_id = tape_cur_id + (1:1:len_v)';

extendtape(T);
h = revADm(v,new_id);

end

%Help function to GAMMA
function y = gamder(x) %used for gamma
if(x<0)
    y = 1./(x.^2).*(x.*gamder(x+1) - gamma(x+1));
else
    y = psi(x)*gamma(x);
end
end

%HELP FUNCTION to "ALL" operations
function [x_type,isAD,n_row,n_col,x_val] = var_type(x)
%VTYPE Return the type of a variable
%   List of posible types:
%   Matrix
%   ColumnVector
%   RowVector
%   Scalar

%% CODE
if (isa(x,'revADm'))
    isAD = true;
    len_x_row = size(x.val,1);
    len_x_col = size(x.val,2);
else
    isAD = false;
    len_x_row = size(x,1);
    len_x_col = size(x,2);
end

if(isAD)
    [n_row,n_col] = size(x.val);
    x_val = x.val;
else
    [n_row,n_col] = size(x);
    x_val = x;
end

if(len_x_row>1 && len_x_col>1)
    x_type = matrixType.Matrix;
elseif(len_x_row>1 && 1==len_x_col)
    x_type = matrixType.ColumnVector;
elseif(1==len_x_row && len_x_col>1)
    x_type = matrixType.RowVector;
elseif(1==len_x_row && 1==len_x_col)
    x_type = matrixType.Scalar;
else
    %ERROR
end
end
function [new_id] = one_arg_ins(old_id,der)
%ONE_ARG_INS takes a objects id and derivative and return the new id
%The IN-parameters:
%OLD_ID     :   is the .id of the object
%DER        :   is the an objects derivative
%The OUT-parameter:
%NEW_id     :   is the .id for the new object


global tape;
global tape_cur_id;

n_el  = sum(sum(old_id>0));     %Number of new elements, i.e. rows in the tape

T     = sparse(n_el,tape_cur_id); % allocate the tape to send

row_id = (1:1:n_el)';           %Create the new rows
col_id_t = old_id(:);           %Create the corresponding columns

col_id   = col_id_t(col_id_t>0);%Only use the ids that is interesting >0

v = tape_cur_id+1:1:tape_cur_id+n_el;

new_id = zeros(numel(old_id),1);                    %Allocate the new_id
new_id(col_id_t>0) = v;                             %Assign the values

new_id  = reshape(new_id,size(der,1),size(der,2));  %reshape the ids to match the output

idx     = sub2ind(size(T),row_id,col_id);           %create a vector of linear index
dert    = der(:);                                   %The values to assign the tape with from the linear index
T(idx)  = dert(col_id_t>0);

extendtape(T);
end
function h = two_arg_interior(u,v,h_val,T1_val,T2_val)
global tape_cur_id;
[utype,is_u,~,~,~] = var_type(u);
[vtype,is_v,~,~,~] = var_type(v);

%CASES decomposition
if(is_u)
    if(matrixType.Scalar==utype)
        if(is_v)
            if(matrixType.Scalar==vtype)
                %CASE 1 - u is revADm and scalar, v is revADm and scalar
                if(max(u.id>0,v.id>0))                                     %hade gått att skriva n=max(u.id>0,v.id>0) men då klagar MATLAB på att det är en bool
                    T = sparse(1,tape_cur_id);
                    if(u.id>0)
                        T(u.id) = T1_val;
                    end
                    if(v.id>0)
                        T(v.id) = T(v.id) + T2_val;
                    end
                    
                    extendtape(T);
                    h = revADm(h_val,tape_cur_id);
                    
                else
                    h = revADm(h_val,0);
                end
            else
                %CASE 2 - u is revADm and scalar, i is revADm but NOT scalar
                if(u.id>0 && max(max(v.id>0)))
                    T = sparse(numel(v.id),tape_cur_id);
                    T(:,u.id) = T1_val(:);
                    
                    v_idt  = v.id;
                    v_id   = v_idt(v_idt>0);
                    
                    n = nnz(v.id);
                    calc_mat = reshape(1:1:n,size(v.id));
                    row_id_t = calc_mat(v.id>0);
                    idx    = sub2ind(size(T),row_id_t,v_id);
                    T2_val_new = ones(size(v.id)).*T2_val;
                    T(idx) = T(idx) + T2_val_new(v.id>0);
                    
                    new_id = row_id_t;
                    new_id(new_id>0)=new_id(new_id>0)+tape_cur_id;
                    
                    extendtape(T);
                    h      = revADm(h_val,reshape(new_id,size(v.id)));
                    
                elseif(u.id>0)
                    T = sparse(numel(v.id),tape_cur_id);
                    T(:,u.id) = T2_val(:);
                    extendtape(T);
                    
                    n = numel(v.id);
                    new_id = zeros(size(v.id));
                    new_id = reshape(tape_cur_id+1:1:tape_cur_id+n,size(new_id));
                    
                    h      = revADm(h_val,new_id);
                elseif(max(v.id>0))
                    v_idt  = v.id;
                    v_id   = v_idt(v_idt>0);
                    n = nnz(v.id);
                    T      = sparse(n,tape_cur_id);
                    idx    = sub2ind(size(T),(1:1:n)',v_id);
                    T(idx) = T1_val;
                    extendtape(T);
                    
                    new_id_vec     = (tape_cur_id+1:1:tape_cur_id+n);
                    new_id = zeros(size(v.id));
                    new_id(v.id>0) = new_id_vec;
                    
                    h = revADm(h_val,new_id);
                    
                else
                    h = revADm(h_val,sparse(size(v.id,1),size(v.id,2)));
                end
            end
        else
            %v is not RevADm
            if(matrixType.Scalar==vtype)
                %CASE 3 - u is revADm and scalar, v is NOT revADm and scalar
                if(u.id>0)
                    T = sparse(1,tape_cur_id);
                    
                    T(1,u.id) = T1_val;
                    extendtape(T);
                    h = revADm(h_val,tape_cur_id);
                else
                    h = revADm(h_val,0);
                end
            else
                %CASE 4 - v is not scalar and not RevADm and u is RevADm and Scalar
                if(u.id>0)
                    n = numel(v);
                    T = sparse(n,tape_cur_id);
                    T(:,u.id) = T1_val(:);
                    
                    new_id = zeros(size(v));
                    new_id = reshape(tape_cur_id+1:1:tape_cur_id+n,size(new_id));
                    
                    extendtape(T);
                    h      = revADm(h_val,new_id);
                else
                    h      = revADm(h_val,zeros(size(v.id)));
                end
            end
        end
    else
        %u is revADm but NOT scalar,
        if(is_v)
            if(matrixType.Scalar==vtype)
                %CASE 5 - u is revADm but NOT scalar, v is RevADm and scalar (compare CASE 2)
                if(v.id>0 && max(max(u.id>0)))
                    T = sparse(numel(u.id),tape_cur_id);
                    T(:,v.id) = T2_val(:);
                    
                    u_idt  = u.id;
                    u_id   = u_idt(u_idt>0);
                    
                    n = nnz(u.id);
                    calc_mat = reshape(1:1:n,size(u.id));
                    row_id_t = calc_mat(u.id>0);
                    idx    = sub2ind(size(T),row_id_t,u_id);
                    T1_val_new = ones(size(u.id)).*T1_val;
                    T(idx) = T(idx) + T1_val_new(u.id>0);
                    
                    new_id = row_id_t;
                    new_id(new_id>0)=new_id(new_id>0)+tape_cur_id;
                    
                    extendtape(T);
                    h      = revADm(h_val,reshape(new_id,size(u.id)));
                elseif(u.id>0)
                    T = sparse(numel(u.id),tape_cur_id);
                    T(:,v.id) = T1_val(:);
                    extendtape(T);
                    
                    n = numel(u.id);
                    new_id = zeros(size(u.id));
                    new_id = reshape(tape_cur_id+1:1:tape_cur_id+n,size(new_id));
                    
                    h      = revADm(h_val,new_id);
                elseif(max(u.id>0))
                    u_idt  = u.id;
                    u_id   = u_idt(u_idt>0);
                    n = nnz(u.id);
                    T      = sparse(n,tape_cur_id);
                    idx    = sub2ind(size(T),(1:1:n)',u_id);
                    T(idx) = T2_val;
                    extendtape(T);
                    
                    new_id_vec     = (tape_cur_id+1:1:tape_cur_id+n);
                    new_id = zeros(size(u.id));
                    new_id(u.id>0) = new_id_vec;
                    
                    h = revADm(h_val,new_id);
                    
                else
                    h = revADm(h_val,sparse(size(u.id,1),size(u.id,2)));
                    
                end
            else
                %CASE 6 - u is RevADm and NOT scalar, v is RevADm and NOT scalar
                if(max(max(u.id>0)) && max(max(v.id>0)))
                    sum_id = u.id+v.id;
                    tot_id = sum_id>0;
                    n = nnz(tot_id);
                    T = sparse(n,tape_cur_id);
                    
                    new_id = zeros(size(tot_id));
                    new_id(tot_id) = 1:1:n;
                    
                    row_id = new_id;
                    
                    new_id(new_id>0) = new_id(new_id>0) + tape_cur_id;
                    
                    u_idt  = u.id;
                    u_id   = u_idt(u_idt>0);
                    
                    v_idt  = v.id;
                    v_id   = v_idt(v_idt>0);
                    
                    %Contribution from u
                    idx1    = sub2ind(size(T),row_id(u.id>0),u_id);
                    to_add1 = T1_val(u.id>0);                        %To get as a column vector
                    T(idx1(:)) = T(idx1(:)) + to_add1(:);
                    %Contribution from v
                    idx2    = sub2ind(size(T),row_id(v.id>0),v_id);
                    to_add2 = T2_val(v.id>0);                        %To get as a column vector
                    T(idx2(:)) = T(idx2(:)) + to_add2(:);
                    
                    extendtape(T);
                    
                    h = revADm(h_val,new_id);
                elseif(max(max(u.id>0)))
                    sum_id = u.id;
                    tot_id = sum_id>0;
                    n = nnz(tot_id);
                    T = sparse(n,tape_cur_id);
                    
                    new_id = zeros(size(tot_id));
                    new_id(tot_id) = 1:1:n;
                    row_id = new_id;
                    
                    new_id(new_id>0) = new_id(new_id>0) + tape_cur_id;
                    
                    u_idt  = u.id;
                    u_id   = u_idt(u_idt>0);
                    idx1    = sub2ind(size(T),row_id(u.id>0),u_id);
                    to_add1 = T1_val(u.id>0);                        %To get as a column vector
                    T(idx1(:)) = T(idx1(:)) + to_add1(:);
                    
                    extendtape(T);
                    h = revADm(h_val,new_id);
                elseif(max(max(v.id>0)))
                    sum_id = v.id;
                    tot_id = sum_id>0;
                    n = nnz(tot_id);
                    T = sparse(n,tape_cur_id);
                    
                    calc_mat = reshape(1:1:numel(v.id),size(v.id));
                    row_id_t = calc_mat(tot_id>0);
                    new_id = row_id_t;
                    new_id(new_id>0)=new_id(new_id>0)+tape_cur_id;
                    new_id = reshape(new_id,size(v.id));
                    
                    v_idt  = v.id;
                    v_id   = v_idt(v_idt>0);
                    
                    %Contribution from v
                    idx2    = sub2ind(size(T),calc_mat(v.id>0),v_id(v.id>0));
                    T(idx2(:)) = T(idx2(:)) + T2_val(:);
                    
                    extendtape(T);
                    h = revADm(h_val,new_id);
                else
                    h = revADm(h_val,zeros(size(u.id)));
                end
            end
        else
            %CASE 7 - u is RevADm NOT scalar, v is not RevADm (scalar or not scalar)
            if(max(max(u.id))>0)
                sum_id = u.id;
                tot_id = sum_id>0;
                n = nnz(tot_id);
                T = sparse(n,tape_cur_id);
                
                calc_mat = reshape(1:1:numel(u.id),size(u.id));
                row_id_t = calc_mat(tot_id>0);
                new_id = row_id_t;
                new_id(new_id>0)=new_id(new_id>0)+tape_cur_id;
                new_id = reshape(new_id,size(u.id));
                
                u_idt  = u.id;
                u_id   = u_idt(u_idt>0);
                
                %Contribution from u
                idx1    = sub2ind(size(T),calc_mat(u.id>0),u_id(u.id>0));
                T(idx1(:)) = T(idx1(:)) + T1_val(:);
                
                extendtape(T);
                h = revADm(h_val,new_id);
            else
                h = revADm(h_val,zeros(size(u.id)));
            end
        end
    end
else
    %u is not RevADm
    if(is_v)
        if(matrixType.Scalar==vtype)
            if(matrixType.Scalar==utype)
                %CASE 8 - u is NOT RevADm but scalar, v is RevADm and scalar
                if(v.id>0)
                    T = sparse(1,tape_cur_id);
                    T(v.id) = T2_val;
                    extendtape(T);
                    
                    h = revADm(h_val,tape_cur_id);
                else
                    h = revADm(h_val,0);
                end
            else
                %CASE 9 - u is NOT RevADm and NOT scalar, v is RevADm and scalar
                if(v.id>0)
                    n = numel(u);
                    T = sparse(n,tape_cur_id);
                    T(:,v.id) = T2_val(:);
                    
                    new_id = zeros(size(u));
                    new_id = reshape(tape_cur_id+1:1:tape_cur_id+n,size(new_id));
                    
                    extendtape(T);
                    h      = revADm(h_val,new_id);
                else
                    h      = revADm(h_val,zeros(size(u)));
                end
            end
        else
            if(max(max(v.id))>0)
                %CASE 10 - u is not RevADm (scalar or not scalar), v is RevADm and NOT scalar
                sum_id = v.id;
                tot_id = sum_id>0;
                n = nnz(tot_id);
                T = sparse(n,tape_cur_id);
                
                v_idt  = v.id;
                v_id   = v_idt(v_idt>0);
                
                new_idt = zeros(size(tot_id));
                new_idt(tot_id) = 1:1:n;
                
                row_id = new_idt;
                
                new_idt = new_idt + tape_cur_id;
                idx2t    = sub2ind(size(T),row_id(v.id>0),v_id);
                to_add2 = T2_val(v.id>0);                        %To get as a column vector
                Tt = T;
                Tt(idx2t(:)) = Tt(idx2t(:)) + to_add2(:);
                
                extendtape(Tt);
                h = revADm(h_val,new_idt);
            else
                %CASE 11 - u is not RevADm (scalar or not scalar), v is NOT RevADm and NOT scalar
                h = revADm(h_val,zeros(size(v.id)));
            end
        end
    else
        error('USER-ERROR in two_arg_interior in AD class');
    end
end
end