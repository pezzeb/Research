classdef forADm
    properties
        val %Function value
        der %derivative value or gradient vector
        %lägga till andra ordningens genom att skapa en ny property der2
        %och sedan överlagra alla andra operationer också
    end
    
    methods
        function obj = forADm(a,b)
            %VALDER class constructor; only the bottom case is needed.
            if nargin == 0 %never intended for use.
                obj.val = [];
                obj.der = [];
            elseif nargin == 1 %c=forADm(a) for constant w/ derivative 0
                obj.val = a;
                obj.der = 0;
            else
                obj.val = a; %given function value
                obj.der = b; %given derivative value or gradient vector
            end
        end
        function vec = double(obj)
            %VALDER/DOUBLE Convert forADm object to vector of doubles.
            vec = [ obj.val, obj.der ];
        end
        function h = plus(u,v)
            %VALDER/MPLUS overloads + with at least one forADm
            if ~isa(u,'forADm') %u is a scalar
                h = forADm(u + v.val, v.der);
            elseif ~isa(v,'forADm') %v is a scalar
                h = forADm(v + u.val, u.der);
            else
                h = forADm(u.val + v.val, u.der + v.der);
            end
        end
        function h = minus(u,v)
            %VALDER/MPLUS overloads + with at least one forADm
            if ~isa(u,'forADm') %u is a scalar
                h = forADm(u - v.val, v.der);
            elseif ~isa(v,'forADm') %v is a scalar
                h = forADm(v - u.val, u.der);
            else
                h = forADm(u.val - v.val, u.der - v.der);
            end
        end
        function h = uminus(u)
            %VALDER/UMINUS overloads + with at least one forADm
            h = forADm(-u.val,-u.der);
        end
        function h = mtimes(u,v)
            %VALDER/MTIMES overloads * with at least one forADm
            if ~isa(u,'forADm') %u is a scalar
                if size(u,1)>1 && size(u,2)>1   %EGEN IF
                    
                end
                h = forADm(u*v.val, u*v.der);
            elseif ~isa(v,'forADm') %v is a scalar
                h = forADm(v*u.val, v*u.der);
            else
                h = forADm(u.val*v.val, u.der*v.val + u.val*v.der);
            end
        end
        function h = times(u,v)
            %VALDER/MTIMES overloads * with at least one forADm
            if ~isa(u,'forADm') %u is a scalar
                h = forADm(u.*v.val, u.*v.der);
            elseif ~isa(v,'forADm') %v is a scalar
                h = forADm(v.*u.val, v.*u.der);
            else
                h = forADm(u.val.*v.val, u.der.*v.val + u.val.*v.der);
            end
        end
        function h = mrdivide(u,v)
            %VALDER/MTIMES overloads / with at least one forADm
            %Täljare/nämnare
            if ~isa(u,'forADm') %u is a scalar
                h = forADm(u/v.val, -u/(v.der)^2);
            elseif ~isa(v,'forADm') %v is a scalar
                h = forADm(u.val/v, u.der/v);
            else
                h = forADm(u.val/v.val, (u.der - (u.val/v.val)*v.der)/v.val);
            end
        end
        function h = power(u,v)
            %VALDER/POWER overloads .^ with at least one forADm
            if ~isa(u,'forADm') %u is a scalar
                %                 h = forADm(u^v.val, u^v.val*log(u)*v.der); TODO
            elseif ~isa(v,'forADm') %v is a scalar
                h = forADm(u.val.^v, v*diag(u.val.^(v-1))*u.der);   %RÄTT?
            else
                h = exp(v*log(u)); %call overloaded log, * and exp
            end
        end
        function h = mpower(u,v)
            %VALDER/MPOWER overloads ^ with at least one forADm
            if ~isa(u,'forADm') %u is a scalar
                h = forADm(u^v.val, u^v.val*log(u)*v.der);
            elseif ~isa(v,'forADm') %v is a scalar
                h = forADm(u.val^v, v*u.val^(v-1)*u.der);
            else
                h = exp(v*log(u)); %call overloaded log, * and exp
            end
        end
        function h = exp(u)
            h = forADm(exp(u.val), exp(u.val)*u.der);
        end
        function h = log(u)
            h = forADm(log(u.val), 1/u.val*u.der);
        end
        function h = sqrt(u)
            h = forADm(sqrt(u.val), (0.5/sqrt(u.val))*u.der);
        end
        function h = sin(u)
            h = forADm(sin(u.val), cos(u.val)*u.der);
        end
        function h = cos(u)
            h = forADm(cos(u.val), -sin(u.val)*u.der);
        end
        function h = tan(u)
            h = forADm(sin(u.val)/cos(u.val), sec(u.val)^2*u.der);
        end
        
        function h = sum(u,dim)
            if nargin==1
                dim = 1;
            end
            h = forADm(sum(u.val,dim), sum(u.der,dim));
        end
        function h = normpdfSPEED(x,mu,sigma)
            
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
                h_val = normpdf(x,mu.val,sigma.val);
                
                x_mu = (x-mu.val);
                h_der = h_val/sigma.val.*(-sigma.der+x_mu/(sigma.val^2).*(mu.der*sigma.val + x_mu*sigma.der));
                
                h = forADm(h_val,h_der);
            catch
                error(message('stats:normpdf:InputSizeMismatch'));
            end
        end
        function h = f_epsilon(X,y)
            Y = X*X';
            
            Yy = (Y.val\y.val);
            
            eps_der_r = -Y.val\(Y.der*Yy) + Y.val\y.der;
            
            R = forADm(Yy,eps_der_r);
            
            h = X'*R;
        end
        function mc = normal_central_moments(sigma,n_m)
            global tape_cur_id;
            switch n_m
                case 1
                    mc = forADm(0,0);
                case 2
                    mc = forADm([0;sigma.val^2],[0;2*sigma.val*sigma.der]);
                case 3
                    mc = forADm([0;sigma.val^2;0],[0;2*sigma.val*sigma.der;0]);
                case 4
                    mc = forADm([0;sigma.val^2;0;3.*sigma.val^4],[0;2*sigma.val*sigma.der;0;12.*sigma.val^3*sigma.der]);
                case 5
                    mc = forADm([0;sigma.val^2;0;3.*sigma.val^4;0],[0;2*sigma.val*sigma.der;0;12.*sigma.val^3;0]);
                case 6
                    mc = forADm([0;sigma.val^2;0;3.*sigma.val^4;0;15.*sigma.val^6],[0;2*sigma.val*sigma.der;0;12.*sigma.val^3;0;90.*sigma.val^6*sigma.der]);
                    %                 case 3
                    %                     new_id  = [0;tape_cur_id+1;0];
                    %
                    %                     T           = sparse(1,tape_cur_id);
                    %                     T(sigma.id) = 2*sigma.val;
                    %                     extendtape(T);
                    %
                    %                     mc = forADm([0;sigma.val^2;0],new_id);
                    %                 case 4
                    %                     new_id  = [0;tape_cur_id+1;0;tape_cur_id+2];
                    %
                    %                     T           = sparse(2,tape_cur_id);
                    %                     T(:,sigma.id) = [2*sigma.val;12*sigma.val^3];
                    %                     extendtape(T);
                    %
                    %                     mc = forADm([0;sigma.val^2;0;3.*sigma.val^4],new_id);
                    %                 case 5
                    %                     new_id  = [0;tape_cur_id+1;0;tape_cur_id+2;0];
                    %
                    %                     T           = sparse(2,tape_cur_id);
                    %                     T(:,sigma.id) = [2*sigma.val;12*sigma.val^3];
                    %                     extendtape(T);
                    %
                    %                     mc = forADm([0;sigma.val^2;0;3.*sigma.val^4;0],new_id);
                    %                 case 6
                    %                     new_id  = [0;tape_cur_id+1;0;tape_cur_id+2;0;tape_cur_id+3];
                    %
                    %                     T           = sparse(3,tape_cur_id);
                    %                     T(:,sigma.id) = [2*sigma.val;12*sigma.val^3;90*sigma.val^5];
                    %                     extendtape(T);
                    %
                    %                     mc = forADm([0;sigma.val^2;0;3.*sigma.val^4;0;15.*sigma.val^6],new_id);
                otherwise
                    error('normal_central_moments error forADm');
            end
        end
        
        function [r,c] = size(u,dim)
            if(exist('dim','var'))
                r = size(u.val,dim);
            else
                [r,c] = size(u.val);
            end
        end
        function h = ctranspose(u)
            h = forADm(u.val',u.der');
        end
        function h = transpose(u)
            h = forADm(u.val.',u.der.');
        end
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
                case '()'   
                    if(length(s.subs)==1)
                        h = forADm(u.val(switch_subs{1}),u.der(switch_subs{1})); %Colon or linearindex
                    else
                        h = forADm(u.val(switch_subs{1,1},switch_subs{1,2}),u.der(switch_subs{1,1},switch_subs{1,2})); %all other for cases of the type (.,.)
                    end     
                case'.'
                    switch s.subs
                        case 'val'
                            h=u.val;
                        case 'der'
                            h=u.der;
                        otherwise
                            error('no such property');
                    end
                case'.()'   %obj.data()
                    switch_subs_val_der = switch_subs{1};
                    switch_subs_rest   = [switch_subs{2:end}];
                    
                    %Determine if val or id (or error)
                    if(strcmp(switch_subs_val_der,'val'))
                        val_id = u.val;
                    elseif(strcmp(switch_subs_val_der,'der'))
                        val_id = u.der;
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
                        
                        h_val = u.val;
                        h_der = u.der;
                        
                        h_val(switch_subs{1}) = b.val;
                        h_der(switch_subs{1}) = b.der;
                        h = forADm(h_val,h_der);
                    else
                            valt = u.val;
                            valt(switch_subs{1},switch_subs{2}) = b.val;
                            
                            dert  = u.der;
                            dert(switch_subs{1},switch_subs{2}) = b.der;
                            
                            h = forADm(valt,dert);
                    end
                otherwise
                    error('error subsasgn');
            end
        end
        
        
    end
end