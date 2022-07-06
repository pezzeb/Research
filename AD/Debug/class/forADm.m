classdef forADm
    properties
        val %Function value
        der %derivative value or gradient vector
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
        %=======OPERATION OF TWO VARIABLES=========================================
        function h = plus(u,v)
            %VALDER/PLUS overloads + with at least one forADm
            [utype,~,urow,ucol,~] = var_type(u);
            [vtype,~,vrow,vcol,~] = var_type(v);
            
            if ~isa(u,'forADm') %u is NOT forADm. All fall finns inte nedan men de som saknas är inte giltiga och matlab kommer kasta fel ändå för dessa
                if(matrixType.Scalar==utype && matrixType.Scalar==vtype)                 %skalär och skalär
                    h = forADm(u + v.val, v.der);
                elseif(matrixType.Scalar==utype && matrixType.ColumnVector==vtype)       %skalär och kolumn
                    h = forADm(u + v.val, v.der);
                elseif(matrixType.Scalar==utype && matrixType.RowVector==vtype)          %skalär och rad
                    h = forADm(u + v.val, v.der);
                elseif(matrixType.Scalar==utype && matrixType.Matrix==vtype)             %skalär och matris
                    h = forADm(u + v.val, v.der);
                    %**********************************************************
                elseif(matrixType.ColumnVector==utype && matrixType.Scalar==vtype)       %kolumn och skalär
                    h = forADm(u + v.val, repmat(v.der,urow,1));
                elseif(matrixType.ColumnVector==utype && matrixType.ColumnVector==vtype) %kolumn och kolumn
                    h = forADm(u + v.val, v.der);
                    %**********************************************************
                elseif(matrixType.RowVector==utype && matrixType.Scalar==vtype)           %rad och skalär
                    h = (u.'+v).';
                elseif(matrixType.RowVector==utype && matrixType.RowVector==vtype)        %rad och rad
                    h = (u.'+v.').';
                elseif(matrixType.RowVector==utype && matrixType.ColumnVector==vtype)        %rad och kolumn
                    Der = [];
                    for ijk=1:ucol
                        Der = [Der;v.der];
                    end
                    h = forADm(u + v.val, Der );
                elseif(matrixType.RowVector==utype && matrixType.Matrix==vtype)        %rad och matrix
                    h = forADm(u+v.val,v.der);
                    %**********************************************************
                elseif(matrixType.Matrix==utype && matrixType.Matrix==vtype)              %matris och matris
                    h = forADm(u + v.val, v.der);
                elseif(matrixType.Matrix==utype && matrixType.ColumnVector==vtype)        %matris och kolumn
                    Der = [];
                    for ijk=1:urow
                        Der = [Der;v.der];
                    end
                    h = forADm(u + v.val, Der);
                    
                else
                    error(message('forADm:plus:Not implemented and perhaps not even defined'));
                end
            elseif ~isa(v,'forADm') %v is NOT forADm - återanvänder kod ovan
                h = v + u;
            else                    % Båda två är forADm 
                %Anledningen till att jag gör så många fall uppdelningar är 
                %för att jag vill få fel för saker som inte är ok i linalg 
                %ex. är rad+kolumn ok i matlab men ej i linalg och ej här.
                if(matrixType.Scalar==utype && matrixType.Scalar==vtype)                  %skalär och skalär
                    h = forADm(u.val + v.val, u.der + v.der);                             
                elseif(matrixType.Scalar==utype && matrixType.ColumnVector==vtype)        %skalär och kolumn
                    h = forADm(u.val + v.val, u.der + v.der);
                elseif(matrixType.Scalar==utype && matrixType.RowVector==vtype)           %skalär och rad
                    h = (u+v.').';
                elseif(matrixType.Scalar==utype && matrixType.Matrix==vtype)              %skalär och matris
                    h = forADm(u.val + v.val, u.der + v.der);
                    %**********************************************************
                elseif(matrixType.ColumnVector==utype && matrixType.Scalar==vtype)        %kolumn och skalär
                    h = v + u;
                elseif(matrixType.ColumnVector==utype && matrixType.ColumnVector==vtype)  %kolumn och kolumn
                    h = forADm(u.val + v.val, u.der + v.der);
                elseif(matrixType.ColumnVector==utype && matrixType.Matrix==vtype)        %kolumn och matris
                    
                    Der = [];
                    vd = v.der;
                    m  = urow;
                    for ijk=1:vcol
                        Der = [Der; vd((ijk-1)*m+1:(ijk*m),:)];
                    end
                    h = forADm(u.val+v.val,Der);
                elseif(matrixType.Matrix==utype && matrixType.ColumnVector==vtype)        %matris och kolumn
                    h = v+u;
                    %**********************************************************
                elseif(matrixType.RowVector==utype && matrixType.Scalar==vtype)           %rad och skalär
                    h = v + u;
                elseif(matrixType.RowVector==utype && matrixType.RowVector==vtype)        %rad och rad
                    h = forADm(u.val + v.val, u.der + v.der);
                    %**********************************************************
                elseif(matrixType.Matrix==utype && matrixType.Matrix==vtype)              %matris och matris
                    h = forADm(u.val + v.val, u.der + v.der);
                else
                    error(message('forADm:plus:Not implemented and perhaps not even defined'));
                end
            end
        end
        function h = minus(u,v)
            %VALDER/MINUS overloads - with at least one forADm
            h = -v+u;
        end
        function h = uminus(u)
            %VALDER/UMINUS overloads - with one forADm
            h = forADm(-u.val,-u.der);
        end
        function h = mtimes(u,v)
            %VALDER/MTIMES overloads * with at least one forADm
            [utype,~,urow,~   ,~] = var_type(u);
            [vtype,~,vrow,vcol,~] = var_type(v);
            if ~isa(u,'forADm')
                if(matrixType.Scalar==utype && matrixType.Scalar==vtype)             %u är skalär och v är skalär
                    h = forADm(u*v.val, u*v.der);
                elseif(matrixType.Scalar==utype && matrixType.RowVector==vtype)      %u är skalär och v är rad
                    h = forADm(u*v.val, u*v.der);
                elseif(matrixType.Scalar==utype && matrixType.ColumnVector==vtype)   %u är skalär och v är kolumn
                    h = (u*(v.')).';
                elseif(matrixType.Scalar==utype && matrixType.Matrix==vtype)         %u är skalär och v är matris
                    vt = forADm(v.val(:),v.der);
                    ht = u*vt;
                    h = forADm(u*v.val,ht.der);
                    
                elseif(matrixType.ColumnVector==utype && matrixType.Scalar==vtype)   %u är kolumn och v är skalär
%                     h = forADm(u*v.val,repmat(u,1,vcol).*repmat(v.der,urow,1)); % not the best testing
                    h = forADm(u*v.val,u.*v.der);
                elseif(matrixType.RowVector==utype && matrixType.Scalar==vtype)      %u är rad    och v är skalär
                    h = (u.'*v).';
                elseif(matrixType.Matrix==utype && matrixType.Scalar==vtype)         %u är matris och v är skalär
                    h = forADm(u*v.val,u(:).*v.der);
                    
                elseif(matrixType.RowVector==utype && matrixType.ColumnVector==vtype)%u är rad och v är kolumn
                    h = forADm(u*v.val, u*v.der);
                elseif(matrixType.ColumnVector==utype && matrixType.RowVector==vtype)%u är kolumn och v är rad                    
                    Der = [];
                    for kl=1:vcol
                        for ij=1:urow
                            dt  = (u(ij)*(v.der(:,kl).'));
                            Der = [Der;dt];
                        end
                    end
                    h = forADm(u*v.val,Der);
                    
                elseif(matrixType.RowVector==utype && matrixType.Matrix==vtype)      %u är rad och v är matris
                    Der = [];
                    m = vrow;
                    for ijk=1:vcol
                        vt = forADm(v.val(:,ijk),v.der((1+(ijk-1)*m):ijk*m,:));
                        vtt=u*vt;
                        Der = [Der,vtt.der.'];
                    end
                    h = forADm(u*v.val,Der);
                    
                elseif(matrixType.Matrix==utype && matrixType.ColumnVector==vtype)   %u är matris och v är kolumn
                    h = forADm(u*v.val, u*v.der);
                elseif(matrixType.Matrix==utype && matrixType.Matrix==vtype)         %u är matris och v är matris
                    Der = [];
                    m = vrow;
                    for ijk=1:vcol
                        vt = forADm(v.val(:,ijk),v.der((1+(ijk-1)*m):ijk*m,:));
                        vtt=u*vt;
                        Der = [Der;vtt.der];
                    end
                    h = forADm(u*v.val,Der);
                else
                    error(message('forADm:mtimes:Not implemented and perhaps not even defined'));
                end
            elseif ~isa(v,'forADm')
                h = (v.'*u.').';
            else %Båda är forADm
                if(matrixType.Scalar==utype && matrixType.Scalar==vtype)             %u är skalär och v är skalär
                    h = forADm(u.val*v.val, u.val*v.der+v.val*u.der);
                elseif(matrixType.Scalar==utype && matrixType.RowVector==vtype)      %u är skalär och v är rad
                    Der = [];
%                     for ijk=1:vcol
%                         Der = [Der, u.val*v.der(:,ijk)+v.val(1,ijk)*(u.der.') ];
%                     end
%                     hnew = forADm(u.val*v.val, Der); %EJ OK!!!!!!!!!
                    h = forADm(u.val*v.val, u.der.'*v.val + u.val*v.der);
                elseif(matrixType.Scalar==utype && matrixType.ColumnVector==vtype)   %u är skalär och v är kolumn
                    h = (u*(v.')).';
                elseif(matrixType.Scalar==utype && matrixType.Matrix==vtype)         %u är skalär och v är matris
%                     warning('Under evaluation')
                    h = forADm(u.val*v.val,u.val*v.der+v.val(:).*u.der);
%                     Der = [];
%                     m   = vrow;
%                     for ijk=1:vcol
%                         Der = [Der;u.der*v.val(:,ijk) + u.val*v.der((ijk-1)*m+1:ijk*m,:)];
%                     end
%                     h = forADm(u.val*v.val,Der);
                    
                elseif(matrixType.ColumnVector==utype && matrixType.Scalar==vtype)   %u är kolumn och v är skalär
                    h = v.*u;
                elseif(matrixType.RowVector==utype && matrixType.Scalar==vtype)      %u är rad    och v är skalär
                    h = (u.'*v).';
                elseif(matrixType.ColumnVector==utype && matrixType.RowVector==vtype)%u är kolumn    och v är rad
                    Der = [];
                    for ijk=1:vcol
                        for lmn=1:urow
                            Der = [Der;u.val(lmn)*(v.der(:,ijk)).' + v.val(ijk)*u.der(lmn,:)];
                        end
                    end
                    
                    h = forADm(u.val*v.val,Der);
                    
                elseif(matrixType.ColumnVector==utype && matrixType.Scalar==vtype)   %u är kolumn och v är skalär
                    
                elseif(matrixType.RowVector==utype && matrixType.ColumnVector==vtype)%u är rad och v är kolumn
                    h = forADm(u.val*v.val, u.val*v.der+v.val.'*u.der.');
                elseif(matrixType.RowVector==utype && matrixType.Matrix==vtype)      %u är rad och v är matris
                    
                    Der = [];
                    m = vrow;
                    for ijk=1:vcol
                       vt = forADm(v.val(:,ijk),v.der((1+(ijk-1)*m):ijk*m,:));
                       vtt=u*vt;
                       Der = [Der,vtt.der.'];
                    end
                    h = forADm(u.val*v.val,Der);
                    
                elseif(matrixType.Matrix==utype && matrixType.ColumnVector==vtype)   %u är matris och v är kolumn
                    h = (v.'*u.').';
                elseif(matrixType.Matrix==utype && matrixType.Matrix==vtype)         %u är matris och v är matris
                    
                    Der = [];
                    m = vrow;
                    for ijk=1:vcol
                        vt = forADm(v.val(:,ijk),v.der((1+(ijk-1)*m):ijk*m,:));
                        vtt=u*vt;
                        Der = [Der;vtt.der];
                    end
                    h = forADm(u.val*v.val,Der);
                    
                else
                    error(message('forADm:mtimes:Not implemented and perhaps not even defined'));
                end
            end
            
        end
        function h = times(u,v)
            %VALDER/TIMES overloads .* with at least one forADm
            [utype,~,~    ,ucol,~] = var_type(u);
            [vtype,~,vrow ,vcol,~] = var_type(v);
            if ~isa(u,'forADm') %u is not forADm
                if(matrixType.Scalar==utype || matrixType.Scalar==vtype)                   % u är skalär ELLER v är skalär - anvnder MTIMES
                    h = u*v;
                
                elseif(matrixType.ColumnVector==utype && matrixType.ColumnVector==vtype)   % u är kolumn och v är kolumn
                    h = forADm(u.*v.val,u.*v.der);
                elseif(matrixType.RowVector==utype && matrixType.RowVector==vtype)         % u är rad och v är rad
                    h = forADm(u.*v.val,u.*v.der);
                elseif(matrixType.Matrix==utype && matrixType.Matrix==vtype)               % u är matris och v är matris
                    vt = forADm(v.val(:),v.der);
                    ht = u(:).*vt;
                    h = forADm(u.*v.val,ht.der);
                elseif(matrixType.Matrix==utype && matrixType.ColumnVector==vtype)          % u är matris och v är kolumn
                    Der = [];
                    for ijk=1:ucol
                        ut = u(:,ijk);
                        ht = ut.*v;
                        Der = [Der;ht.der];
                    end
                    h = forADm(u.*v.val,Der);
                elseif(matrixType.ColumnVector==utype && matrixType.Matrix==vtype)          % u är kolumn och v är matris    
                    Der = [];
                    m = vrow;
                    for ijk=1:vcol
                       vt = forADm(v.val(:,ijk),v.der( ((ijk-1)*m + 1):((ijk)*m),:));
                       ht = u.*vt;
                       Der = [Der;ht.der];
                    end
                    
                    h = forADm(u.*v.val,Der);
                elseif(matrixType.ColumnVector==utype && matrixType.RowVector==vtype)      % u är kolumn och v är rad
                    Der = [];
                    for ijk=1:vcol
                        Der = [Der;u.*v.der(:,ijk).'];
                    end
                    h = forADm(u.*v.val,Der);
                elseif(matrixType.RowVector==utype && matrixType.ColumnVector==vtype)      % u är rad och v är kolumn
                    Der = [];
                    for ijk=1:ucol
                       Der = [Der;u(ijk)*v.der];
                    end
                    h = forADm(u.*v.val,Der);
                    
                else
                    error(message('forADm:times:Not implemented and perhaps not even defined'));
                end
            elseif ~isa(v,'forADm') %v is not forADm
                h = v.*u;
            else %Båda är forADm
                if(matrixType.Scalar==utype || matrixType.Scalar==vtype)                   %u är skalär ELLER v är skalär
                    h=u*v;
                    
                elseif(matrixType.ColumnVector==utype && matrixType.ColumnVector==vtype)   %u är kolumn och v är kolumn
                    %                     h = forADm(u.val.*v.val,repmat(u.val,1,size(v.der,2)).*v.der + repmat(v.val,1,size(u.der,2)).*u.der );
                    h = forADm(u.val.*v.val,u.val.*v.der + v.val.*u.der );
                    
                elseif(matrixType.RowVector==utype && matrixType.RowVector==vtype)         %u är rad och v är rad
                    h = (u.'.*v.').';
                    
                elseif(matrixType.Matrix==utype && matrixType.Matrix==vtype)               %u är matris och v är matris
                    Der = [];
                    m = vrow;
                    for ijk=1:ucol
                        ut = forADm(u.val(:,ijk),u.der((1+(ijk-1)*m):ijk*m,:));
                        vt = forADm(v.val(:,ijk),v.der((1+(ijk-1)*m):ijk*m,:));
                        uvt = ut.*vt;
                        Der = [Der;uvt.der];
                    end
                    
                    h = forADm(u.val.*v.val,Der);
                
                elseif(matrixType.ColumnVector==utype && matrixType.Matrix==vtype)         %u är kolumn och v är matris
                    Der = [];
                    m = vrow;
                    for ijk=1:vcol
                        vt = forADm(v.val(:,ijk),v.der((1+(ijk-1)*m):ijk*m,:));
                        uvt = u.*vt;
                        Der = [Der;uvt.der];
                    end
                    h = forADm(u.val.*v.val,Der);
                    
                else
                    error(message('forADm:times:Not implemented and perhaps not even defined'));
                end
            end
        end
        function h = mrdivide(u,v)
            %VALDER/MRDIVIDE overloads / with at least one forADm
            %Täljare/nämnare
            [vtype,~,~,~,~] = var_type(v);
            if(matrixType.Scalar==vtype)
                h = u./v;
            else
                error(message('forADm:mrdivide:Not implemented and perhaps not even defined'));
            end
        end
        function h = rdivide(u,v)
            %VALDER/RDIVIDE overloads ./ with at least one forADm
            %Täljare/nämnare
            [utype,~,urow,~   ,~] = var_type(u);
            [vtype,~,vrow,vcol,~] = var_type(v);
            if ~isa(u,'forADm') %u is NOT forADm
                if(matrixType.Scalar==utype && matrixType.Scalar==vtype)                %u är skalär och v är skalär
                    h = forADm(u./v.val, -u./((v.val).^2).*v.der);
                elseif(matrixType.Scalar==utype && matrixType.ColumnVector==vtype)      %u är skalär och v är kolumn
%                     h = forADm(u./v.val, diag(-u./((v.val).^2))*v.der);
                    h = forADm(u./v.val,-u./(v.val.^2).*v.der); 
                elseif(matrixType.Scalar==utype && matrixType.RowVector==vtype)         %u är skalär och v är rad
                    h = (u./(v.')).';
                elseif(matrixType.Scalar==utype && matrixType.Matrix==vtype)            %u är skalär och v är matris
                    vt = forADm(v.val(:),v.der);
                    ht = u./vt;
                    h = forADm(u./v.val,ht.der);
                    %hnew = forADm(u./v.val,-u./(v.val(:).^2).*v.der); %IMPROVMENT? CHECK TODO!
                    
                elseif(matrixType.ColumnVector==utype && matrixType.Scalar==vtype)      % u är kolumn och v är skalär
                    h = forADm(u./v.val,diag(-u./(v.val^2))*repmat(v.der,urow,1));
                    %hnew = forADm(u./v.val,-u/(v.val^2).*v.der);
                elseif(matrixType.ColumnVector==utype && matrixType.ColumnVector==vtype)% u är kolumn och v är kolumn
                    h = forADm(u./v.val,diag(-u./(v.val.^2))*v.der);
                    %hnew = forADm(u./v.val,-u./(v.val.^2).*v.der); %IMPROVMENT? CHECK TODO!
                    
                elseif(matrixType.RowVector==utype && matrixType.Scalar==vtype)         % u är rad och v är skalär
                    h = (u.'./v).';
                elseif(matrixType.RowVector==utype && matrixType.RowVector==vtype)      % u är rad och v är rad
                    h = (u.'./(v.')).';
                    
                elseif(matrixType.Matrix==utype && matrixType.Scalar==vtype)            % u är matris och v är skalär
                    derval = -u./(v.val^2);
                    h = forADm(u./v.val,derval(:).*v.der);
                    % hnew = forADm(u./v.val,-u(:)./(v.val^2).*v.der); %IMPROVMENT? CHECK TODO!
                elseif(matrixType.Matrix==utype && matrixType.Matrix==vtype)            % u är matris och v är matris
                    vt = forADm(v.val(:),v.der);
                    ht = u(:)./vt;
                    h = forADm(u./v.val,ht.der);
                    %hnew = forADm(u./v.val,-u(:)./(v.val(:).^2).*v.der); %IMPROVMENT? CHECK TODO!
                else
                    error(message('forADm:rdivide:Not implemented and perhaps not even defined'));
                end
            elseif ~isa(v,'forADm') %v is not forADm
                 hold = u.*(1./v); %använda denna eller ifsatsen
                if(matrixType.Scalar==utype && matrixType.Scalar==vtype)                %u är skalär och v är skalär
                    h = forADm(u.val/v,u.der/v);
                elseif(matrixType.Scalar==utype && matrixType.ColumnVector==vtype)      %u är skalär och v är kolumn
                    h = forADm(u.val./v,u.der./v);
                elseif(matrixType.Scalar==utype && matrixType.RowVector==vtype)         %u är skalär och v är rad
                    h = (u./(v.')).';
                elseif(matrixType.Scalar==utype && matrixType.Matrix==vtype)            %u är skalär och v är matris
                    h = forADm(u.val./v,u.der./(v(:)));
                    
                elseif(matrixType.ColumnVector==utype && matrixType.Scalar==vtype)      % u är kolumn och v är skalär
                    h = forADm(u.val./v,u.der/v);
                elseif(matrixType.ColumnVector==utype && matrixType.ColumnVector==vtype)% u är kolumn och v är kolumn
                    h = forADm(u.val./v,u.der./v);
                    
                elseif(matrixType.RowVector==utype && matrixType.Scalar==vtype)         % u är rad och v är skalär
                    h = (u.'./v).';
                elseif(matrixType.RowVector==utype && matrixType.RowVector==vtype)      % u är rad och v är rad
                    h = (u.'./(v.')).';
                    
                elseif(matrixType.Matrix==utype && matrixType.Scalar==vtype)            % u är matris och v är skalär
                    h = forADm(u.val/v,u.der/v);
                elseif(matrixType.Matrix==utype && matrixType.Matrix==vtype)            % u är matris och v är matris
                    h = forADm(u.val./v,u.der./v(:));
                else
                    ht = u.*(1./v); %övriga fall detta är inte bra ur precission synpunkt
                    h = forADm(u.val./v,ht.der);
                end
                    
            else
                h = u.*(1./v); %INTE DET SNYGGAST, evnetuellt blir prestandan: snabbhet och precission sämre
            end
        end
        function h = power(u,v)
            %VALDER/POWER overloads .^ with at least one forADm
            %             h = exp(v.*log(u)); %Rent matematiskt kan man göra så här på alla men precissionen blir sämre!!!!!
            [utype,~,urow,~,~] = var_type(u);
            [vtype,~,vrow,~,~] = var_type(v);
            
            if ~isa(u,'forADm') %u is not forADm
                %                 h = forADm(u.^v.val, diag(u.^v.val.*log(u))*v.der);
                if(matrixType.Scalar==utype && matrixType.Scalar==vtype)                 %u är skalär och v är skalär
                    h = forADm(u.^v.val,u.^v.val.*log(u).*v.der);
                elseif(matrixType.Scalar==utype && matrixType.RowVector==vtype)          %u är skalär och v är rad
                    h = (u.^(v.')).';
                elseif(matrixType.Scalar==utype && matrixType.ColumnVector==vtype)       %u är skalär och v är kolumn
                    h = forADm(u.^v.val,diag(u.^v.val.* log(u))*v.der);
                elseif(matrixType.Scalar==utype && matrixType.Matrix==vtype)             %u är skalär och v är matris
                    vt = forADm(v.val(:),v.der);
                    ht = u.^vt;
                    h = forADm(u.^v.val,ht.der);
                    
                elseif(matrixType.ColumnVector==utype && matrixType.Scalar==vtype)       %u är kolumn och v är skalär
                    h = forADm(u.^v.val,diag(u.^v.val.*log(u))*repmat(v.der,urow,1) );
                elseif(matrixType.RowVector==utype && matrixType.Scalar==vtype)          %u är rad    och v är skalär
                    h = (u.'.^v).';
                elseif(matrixType.Matrix==utype && matrixType.Scalar==vtype)             %u är matris och v är skalär
                    ht = u(:).^v;
                    h = forADm(u.^v.val,ht.der);
                    
                elseif(matrixType.ColumnVector==utype && matrixType.ColumnVector==vtype) %u är kolumn och v är kolumn
                    h = forADm(u.^v.val,diag(u.^v.val.* log(u))*v.der);
                elseif(matrixType.RowVector==utype && matrixType.RowVector==vtype)       %u är rad och v är rad
                    h = (u.'.^(v.')).';
                elseif(matrixType.Matrix==utype && matrixType.Matrix==vtype)             %u är matris och v är matris
                    vt = forADm(v.val(:),v.der);
                    ht = u(:).^vt;
                    h = forADm(u.^v.val,ht.der);
                    
                else
                    error(message('forADm:mpower:Not implemented and perhaps not even defined'));
                end
                
            elseif ~isa(v,'forADm') %v is not forADm
                if(matrixType.Scalar==utype && matrixType.Scalar==vtype)                 %u är skalär och v är skalär
                    h = forADm(u.val.^v,v.*u.val.^(v-1).*u.der);
                elseif(matrixType.Scalar==utype && matrixType.RowVector==vtype)          %u är skalär och v är rad
                    h = (u.^(v.'))';
                elseif(matrixType.Scalar==utype && matrixType.ColumnVector==vtype)       %u är skalär och v är kolumn
                    h = forADm(u.val.^v,diag(v.*u.val.^(v-1))* repmat(u.der,vrow,1));
                elseif(matrixType.Scalar==utype && matrixType.Matrix==vtype)             %u är skalär och v är matris
                    ht = u.^v(:);
                    h = forADm(u.val.^v,ht.der);
                    
                elseif(matrixType.ColumnVector==utype && matrixType.Scalar==vtype)       %u är kolumn och v är skalär
                    h = forADm(u.val.^v,v.*u.val.^(v-1).*u.der);
                elseif(matrixType.RowVector==utype && matrixType.Scalar==vtype)          %u är rad    och v är skalär
                    h = (u.'.^v).';
                elseif(matrixType.Matrix==utype && matrixType.Scalar==vtype)             %u är matris och v är skalär
                    ut = forADm(u.val(:),u.der);
                    ht = ut.^v;
                    h = forADm(u.val.^v,ht.der);
                    
                elseif(matrixType.ColumnVector==utype && matrixType.ColumnVector==vtype) %u är kolumn och v är kolumn
                    h = forADm(u.val.^v,diag(v.*u.val.^(v-1))*u.der);
                elseif(matrixType.RowVector==utype && matrixType.RowVector==vtype)       %u är rad och v är rad
                    h = (u.'.^(v.')).';
                elseif(matrixType.Matrix==utype && matrixType.Matrix==vtype)             %u är matris och v är matris
                    ut = forADm(u.val(:),u.der);
                    ht = ut.^v(:);
                    h = forADm(u.val.^v,ht.der);
                    
                else
                    error(message('forADm:mpower:Not implemented and perhaps not even defined'));
                end
            else
                h = exp(v.*log(u)); %call overloaded log, * and exp - Det är rätt matematiskt men det blir kanske sämre precission
            end
        end
        function h = mpower(u,v)
            %VALDER/MPOWER overloads ^ with at least one forADm
            [utype,~,~,~,~] = var_type(u);
            [vtype,~,~,~,~] = var_type(v);
            if(matrixType.Scalar==utype && matrixType.Scalar==vtype)
                h=u.^v;
            elseif(matrixType.Matrix==utype && matrixType.Scalar==vtype)
                if(~isa(v,'forADm'))
                    if(mod(v,1)<1e-10) %slippa numeriska problem
                        ht =  u;
                        for i=1:round(v)-1
                           ht = ht*u;
                        end
                        h=ht;
                    else
                        error(message('forADm:mpower:Not implemented and perhaps not even defined(1)'));
                    end
                else
                    error(message('forADm:mpower:Not implemented and perhaps not even defined(2)'));
                end
            else
                %No other case is interesting
                error(message('forADm:mpower:Not implemented and perhaps not even defined(3)'));
            end
        end
        %Operation of single variable
        function h = exp(u)
            [utype,~,~,~,~] = var_type(u);
            
            if(matrixType.Scalar==utype)
                %skalär
                h = forADm(exp(u.val), exp(u.val)*u.der);
            elseif(matrixType.ColumnVector==utype)
                %kolumnvektor
                h = forADm(exp(u.val), exp(u.val).*u.der);
            elseif(matrixType.RowVector==utype)
                %radvektor
                h = forADm(exp(u.val), exp(u.val).*u.der);
            elseif(matrixType.Matrix==utype)
                %Matris
                h = forADm(exp(u.val),exp(u.val(:)).*u.der);
            end
        end
        function h = log(u)
            [utype,~,~,~,~] = var_type(u);
            if(matrixType.Scalar==utype)
                h = forADm(log(u.val), 1./u.val*u.der);
            elseif(matrixType.ColumnVector==utype)
                h = forADm(log(u.val), 1./u.val.*u.der);
            elseif(matrixType.RowVector==utype)
                h = forADm(log(u.val), u.der*diag(1./u.val));
            elseif(matrixType.Matrix==utype)
                h = forADm(log(u.val),1./u.val(:).*u.der);
            end
        end
        function h = sqrt(u)
            [utype,~,~,~,~] = var_type(u);
            if(matrixType.Scalar==utype)
                h = forADm(sqrt(u.val), 0.5./sqrt(u.val).*u.der);
            elseif(matrixType.ColumnVector==utype)
                h = forADm(sqrt(u.val), 0.5./sqrt(u.val).*u.der);
            elseif(matrixType.RowVector==utype)
                h = forADm(sqrt(u.val), u.der*diag(0.5./sqrt(u.val)));
            elseif(matrixType.Matrix==utype)
                h = forADm(sqrt(u.val), 1./sqrt(u.val(:)).*u.der );
            end
        end
        function h = sin(u)
            [utype,~,~,~,~] = var_type(u);
            
            if(matrixType.Scalar==utype)
                h = forADm(sin(u.val), cos(u.val)*u.der);
            elseif(matrixType.ColumnVector==utype)
                h = forADm(sin(u.val), diag(cos(u.val))*u.der);
            elseif(matrixType.RowVector==utype)
                h = forADm(sin(u.val), u.der*diag(cos(u.val)));
            elseif(matrixType.Matrix==utype)
                h = forADm(sin(u.val), cos(u.val(:)).*u.der);
            end
        end
        function h = cos(u)
            [utype,~,~,~,~] = var_type(u);
            
            if(matrixType.Scalar==utype)
                h = forADm(cos(u.val), -sin(u.val)*u.der);
            elseif(matrixType.ColumnVector==utype)
                h = forADm(cos(u.val), diag(-sin(u.val))*u.der);
            elseif(matrixType.RowVector==utype)
                h = forADm(cos(u.val), u.der*diag(-sin(u.val)));
            elseif(matrixType.Matrix==utype)
                h = forADm(cos(u.val), -sin(u.val(:)).*u.der);
            end
        end
        function h = tan(u)
            %             h = forADm(sin(u.val)/cos(u.val), sec(u.val)^2*u.der);
            h = sin(u)./cos(u);
        end
        function h = conj(z)
            h = forADm(conj(z.val),conj(z.der));
        end
        function h = real(z)
            w = conj(z);
            h = (z + w)/2;
        end
        function h = imag(z)
            h = (z - conj(z))/(2*1i);
        end
        function h = sum(u,dim)
            [utype,~,u_row,~,~] = var_type(u);
            if nargin==1
                dim = 1;
            end
            %Observera skillnadnen mellan if-satserna i if och else-delen
            if(1==dim)
                if(matrixType.Scalar==utype || matrixType.ColumnVector==utype)
                    h = forADm(sum(u.val,dim), sum(u.der,dim));
                elseif(matrixType.RowVector==utype)
                    h = sum(u.',1).';
                elseif(matrixType.Matrix==utype)
                    m = u_row;
                    n = size(u.der,2);
                    nel=numel(u.der);
                    rder = reshape(u.der,m,nel/m);
                    srder=sum(rder,1);
                    rsrder=reshape(srder,m,n);
                    
                    h = forADm(sum(u.val),rsrder.');
                end
            elseif(2==dim)
                if(matrixType.Scalar==utype || matrixType.RowVector==utype)
                    h = forADm(sum(u.val,dim), sum(u.der,dim));
                elseif(matrixType.ColumnVector==utype)
                    h = sum(u.',1).';
                elseif(matrixType.Matrix==utype)
                    h = sum(u.',1).';
                end
            end
        end
        function h = gamma(u)
            fprintf('Är detta verkligen korrekt implemenerat?')
            [utype,~,~,~,~] = var_type(u);
            if(matrixType.Scalar==utype)
                h = forADm(gamma(u.val),gamma(u.val).*psi(u.val).*u.der);
            elseif(matrixType.ColumnVector==utype)
                h = forADm(gamma(u.val),gamma(u.val).*psi(u.val).*u.der);
            elseif(matrixType.RowVector==utype)
                h = forADm(gamma(u.val),gamma(u.val).*psi(u.val).*u.der);
            elseif(matrixType.Matrix==utype)
                h = forADm(gamma(u.val),gamma(u.val(:)).*psi(u.val(:)).*u.der);
            end
        end
        
        function h = normcdf(u,mu,sigma) %Lite "fuskat" eftersom att man måste ange mu och sigma
            [utype,~,~,~,~] = var_type(u);
            if(matrixType.Scalar==utype)
                h = forADm(normcdf(u.val,mu,sigma),normpdf(u.val,mu,sigma).*u.der);
            elseif(matrixType.ColumnVector==utype)
                fprintf('Ej testad')
                h = forADm(normcdf(u.val,mu,sigma), diag(normpdf(u.val,mu,sigma))*u.der);
            elseif(matrixType.RowVector==utype)
                fprintf('Ej implementerat')
                h = forADm(normcdf(u.val,mu,sigma),u.der*diag(normpdf(u.val,mu,sigma)));
            elseif(matrixType.Matrix==utype)
                fprintf('Ej implementerat')
                h = forADm(normcdf(u.val,mu,sigma), normpdf(u.val).*u.der);
            end
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
        function h = fft(u)
            h = forADm(fft(u.val),fft(u.der));
        end
        
        function [r,c] = size(u,dim)
            %Tänk på att denna anropas om man inte har semikolon och resultat i kommandfönstret är forADm
            if(exist('dim','var'))
                r = size(u.val,dim);
            else
                [r,c] = size(u.val);
            end
        end
        function h = ctranspose(u)
            [utype,~,u_row,u_col,~] = var_type(u);
            if(matrixType.Scalar==utype)
                h = forADm(u.val,u.der);
            elseif(matrixType.ColumnVector==utype || matrixType.RowVector==utype)
                h = forADm(u.val',u.der');
            elseif(matrixType.Matrix==utype)
                nel = u_row*u_col; %number of elements
                I = reshape((1:1:nel),u_row,u_col);
                h = forADm(u.val',conj(u.der(I',:))); %conjugate since this is the complex conjugate transpose
            end
        end
        function h = transpose(u)
            [utype,~,u_row,u_col,~] = var_type(u);
            if(matrixType.Scalar==utype)
                h = forADm(u.val,u.der);
            elseif(matrixType.ColumnVector==utype || matrixType.RowVector==utype)
                h = forADm(u.val.',u.der.');
            elseif(matrixType.Matrix==utype)
                nel = u_row*u_col; %number of elements
                I = reshape((1:1:nel),u_row,u_col);
                h = forADm(u.val.',u.der(I',:));
            end
        end
        
        %% OBSERVERA ATT DESSA ENBART GER ETT VÄRDE DÅ DERIVATAN INTE ÄR DEFINERAD
        function hval = round(u)
            hval = round(u.val);
        end
        function hval = floor(u)
            hval = floor(u.val);
        end
        
        %% Logic Operators - (https://se.mathworks.com/help/matlab/matlab_oop/implementing-operators-for-your-class.html?s_tid=gn_loc_drop)
        function hlogi = lt(u,v)
            [uval,vval] = getval(u,v);
            if(uval < vval)
                hlogi = true;
            else
                hlogi = false;
            end
        end
        function hlogi = gt(u,v)
            [uval,vval] = getval(u,v);
            if(uval >  vval)
                hlogi = true;
            else
                hlogi = false;
            end
        end
        function hlogi = le(u,v)
            [uval,vval] = getval(u,v);
            if(uval <= vval)
                hlogi = true;
            else
                hlogi = false;
            end
        end
        function hlogi = ge(u,v)
            [uval,vval] = getval(u,v);
            if(uval >=  vval)
                hlogi = true;
            else
                hlogi = false;
            end
        end
        function hlogi = ne(u,v)
            [uval,vval] = getval(u,v);
            if(uval ~= vval)
                hlogi = true;
            else
                hlogi = false;
            end
        end
        function hlogi = eq(u,v)
            [uval,vval] = getval(u,v);
            if(uval == vval)
                hlogi = true;
            else
                hlogi = false;
            end
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
                        h = forADm(u.val(switch_subs{1,1},switch_subs{1,2}),u.der(switch_subs{1,1},:)); %all other for cases of the type (.,.)
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

%HELP FUNCTIONS
function [x_type,isAD,n_row,n_col,x_val] = var_type(x)
%VTYPE Return the type of a variable
%   List of posible types:
%   Matrix
%   ColumnVector
%   RowVector
%   Scalar

%% CODE
if (isa(x,'forADm'))
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

function [uval,vval] = getval(u,v)

if(~isa(u,'forADm')) % u är inte forADm
    uval = u;
    vval = v.val;
elseif(~isa(v,'forADm')) % u är inte forADm
    uval = u.val;
    vval = v;
else
    uval = u.val;
    vval = v.val;
end

end