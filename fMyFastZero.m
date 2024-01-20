function [b,fval,exitflag,output] = fMyFastZero(FunFcnIn,x,options,varargin)

% Initialization
fcount = 0;
iter = 0;
intervaliter = 0;
exitflag = 1;
procedure = ' ';

% % Set up default options
defaultopt.Display = 'notify';
defaultopt.TolX = eps;
defaultopt.FunValCheck  = 'off';
defaultopt.OutputFcn = [];
defaultopt.PlotFcns = [];

allOptionsDefault = false;

if ~allOptionsDefault
    tol = optimget(options,'TolX',defaultopt,'fast');

else 
    tol = eps;
    trace = 1;
end

if isa(FunFcnIn,'function_handle')
    FunFcn = FunFcnIn;
else
    % Convert to function handle as needed.
    [FunFcn,errStruct] = fcnchk(FunFcnIn,length(varargin));  %#ok<DFCNCHK>
    if ~isempty(errStruct)
        error(message(errStruct.identifier));
    end
end

if ~allfinite(x) % ~all(isfinite(x))
    error('MATLAB:fzero:Arg2NotFinite',...
        getString(message('MATLAB:optimfun:fzero:Arg2NotFinite')));
end

% Interval input
if (numel(x) == 2) 
    a = x(1); savea=a;
    b = x(2); saveb=b;
    % Put first feval in try catch
    fa = localFirstFcnEval(FunFcn,FunFcnIn,a,varargin{:});
    fb = FunFcn(b,varargin{:});
    if ~isfinite(fa) || ~isfinite(fb) || ~isreal(fa) || ~isreal(fb)
        error('MATLAB:fzero:ValuesAtEndPtsComplexOrNotFinite',...
            getString(message('MATLAB:optimfun:fzero:ValuesAtEndPtsComplexOrNotFinite')));
    end
    fcount = fcount + 2;
    savefa = fa; savefb = fb;
    
    if ( fa == 0 )
        b = a;
        fval = fa;
        return
    elseif ( fb == 0)
        % b = b;
        fval = fb;
        return
    elseif (fa > 0) == (fb > 0)
        error('MATLAB:fzero:ValuesAtEndPtsSameSign',...
            getString(message('MATLAB:optimfun:fzero:ValuesAtEndPtsSameSign')));
    end
    
    % Starting guess scalar input
elseif (numel(x) == 1)
    % Put first feval in try catch
    fx = localFirstFcnEval(FunFcn,FunFcnIn,x,varargin{:});
    fcount = fcount + 1;  
    if fx == 0
        b = x;
        fval = fx;
        return
    elseif ~isfinite(fx) || ~isreal(fx)
        error('MATLAB:fzero:ValueAtInitGuessComplexOrNotFinite',...
            getString(message('MATLAB:optimfun:fzero:ValueAtInitGuessComplexOrNotFinite')));
    end
    
    if x ~= 0
        dx = x/50;
    else 
        dx = 1/50;
    end
    
    % Find change of sign.
    twosqrt = sqrt(2); 
    a = x; fa = fx; b = x; fb = fx;

    while (fa > 0) == (fb > 0)
        intervaliter = intervaliter + 1;
        dx = twosqrt*dx;
        a = x - dx;  fa = FunFcn(a,varargin{:});
        fcount = fcount + 1;
%         if ~isfinite(fa) || ~isreal(fa) || ~isfinite(a)
%             [b,fval,exitflag,output] = errHandler(a,fa,intervaliter,iter,fcount,trace,buildOutputStruct);
%             return
%         end

        if (fa > 0) ~= (fb > 0) % check for different sign
            % Before we exit the while loop, print out the latest interval
            break
        end
        
        b = x + dx;  fb = FunFcn(b,varargin{:});
%         if ~isfinite(fb) || ~isreal(fb) || ~isfinite(b)
%             [b,fval,exitflag,output] = errHandler(b,fb,intervaliter,iter,fcount,trace,buildOutputStruct);
%             return
%         end
        fcount = fcount + 1;        

    end % while
    savea = a; savefa = fa; saveb = b; savefb = fb;
else
    error('MATLAB:fzero:LengthArg2',...
        getString(message('MATLAB:optimfun:fzero:LengthArg2')));
end % if (numel(x) == 2)

fc = fb;
procedure = 'initial';

% Main loop, exit from middle of the loop
while fb ~= 0 && a ~= b
    % Insure that b is the best result so far, a is the previous
    % value of b, and c is on the opposite side of the zero from b.
    if (fb > 0) == (fc > 0)
        c = a;  fc = fa;
        d = b - a;  e = d;
    end
    if abs(fc) < abs(fb)
        a = b;    b = c;    c = a;
        fa = fb;  fb = fc;  fc = fa;
    end
    
    % Convergence test and possible exit
    m = 0.5*(c - b);
    toler = 2.0*tol*max(abs(b),1.0);
    if (abs(m) <= toler) || (fb == 0.0) 
        break
    end

    % Choose bisection or interpolation
    if (abs(e) < toler) || (abs(fa) <= abs(fb))
        % Bisection
        d = m;  e = m;
        procedure='bisection';
    else
        % Interpolation
        s = fb/fa;
        if (a == c)
            % Linear interpolation
            p = 2.0*m*s;
            q = 1.0 - s;
        else
            % Inverse quadratic interpolation
            q = fa/fc;
            r = fb/fc;
            p = s*(2.0*m*q*(q - r) - (b - a)*(r - 1.0));
            q = (q - 1.0)*(r - 1.0)*(s - 1.0);
        end
        if p > 0
            q = -q;
        else
            p = -p;
        end
        % Is interpolated point acceptable
        if (2.0*p < 3.0*m*q - abs(toler*q)) && (p < abs(0.5*e*q))
            e = d;  d = p/q;
            procedure='interpolation';
        else
            d = m;  e = m;
            procedure='bisection';
        end
    end % Interpolation
    
    % Next point
    a = b;
    fa = fb;
    if abs(d) > toler
        b = b + d;
    elseif b > c
        b = b - toler;
    else
        b = b + toler;
    end
    fb = FunFcn(b,varargin{:});
    fcount = fcount + 1;
    iter = iter + 1;
end % Main loop

fval = fb; % b is the best value

if abs(fval) > max(abs(savefa),abs(savefb))
    exitflag = -5;
end

%--------------------------------------------------------------------------
function fx = localFirstFcnEval(FunFcn,FunFcnIn,x,varargin)

% Put first feval in try catch
try
    fx = FunFcn(x,varargin{:});
catch ME
    % Additional processing for error handling
    if isa(FunFcn,'function_handle') % function handle
        Ffcnstr = func2str(FunFcn);  % get the name passed in
        Ftype = 'function_handle';
    elseif isa(FunFcn,'inline')
        if isa(FunFcnIn,'inline')
            Ffcnstr = inputname(1);  % name of inline object such as f where f=inline('x*2');
            if isempty(Ffcnstr)  % inline('sin(x)')
                Ffcnstr = formula(FunFcn);  % Grab formula, no argument name
            end
            Ftype = 'inline object';
        else  % not an inline originally (character array expression).
            Ffcnstr = char(FunFcnIn);  % get the character array expression
            Ftype = 'expression';
        end
    else  % Not converted, must be m-file or builtin
        Ffcnstr = char(FunFcnIn);  % get the name passed in
        Ftype = 'function';
    end
    if ~isempty(Ffcnstr)
        error('MATLAB:fzero:InvalidFunctionSupplied',...
            getString(message('MATLAB:optimfun:fzero:InvalidFunctionSupplied',sprintf('%s ==> %s',Ftype,Ffcnstr),ME.message)));
    else
        error('MATLAB:fzero:InvalidFunctionSupplied',...
            getString(message('MATLAB:optimfun:fzero:InvalidFunctionSupplied',Ftype,ME.message)));
    end
    
end