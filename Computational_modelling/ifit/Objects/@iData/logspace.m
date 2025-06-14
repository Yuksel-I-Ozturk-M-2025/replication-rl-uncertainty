function v = logspace(a,b,n)
% v = logspace(a,b,n) : creates a logarithmically spaced vector of objects
%
%   @iData/logspace function to create a logarithmically spaced vector of iData objects
%   This corresponds to a 'morphing' from a to b in n steps.
%
% input:  a: object  (iData)
%         b: object  (iData)
%         n: number of steps, i.e. length of vector. Default is n=10 (integer)
% output: v: vector (iData array)
% ex:     b=logspace(a,b);
%
% Version: Aug. 22, 2017
% See also iData, iData/max, iData/min, iData/colon, iData/linspace

if nargin <= 2
  n = [];
end
if nargin == 1, b=a; end

if isempty(n) || n <=0
  n=10;
end

if isa(a, 'iData') && numel(a) > 1, a=a(1); end
if isa(b, 'iData') && numel(b) > 1, b=b(end); end

% get scalar value (if specified as number of 1x1 object)
if     isempty(a), a=0; 
elseif isempty(b), b=0; 
elseif isscalar(a) & isa(a, 'iData'), a=get(a,'Signal');
elseif isscalar(b) & isa(b, 'iData'), b=get(b,'Signal');
end

% one of (a,b) must be an iData else fallback to the usual linspace
if ~isa(a,'iData') && ~isa(b,'iData')
  v=logspace(a,b,n);
  return
end

% create constant objects from scalar input, using the other object as template
if ~isa(a, 'iData') & isnumeric(a) & ~isempty(b)
  s=mean(a(:))*ones(size(get(b,'Signal'))); a=copyobj(b); set(a,'Signal', s);
elseif ~isa(b, 'iData') & isnumeric(b) & ~isempty(a)
  s=mean(b(:))*ones(size(get(a,'Signal'))); b=copyobj(a); set(b,'Signal', s);
else
  [a,b] = intersect(a,b); % get intersection for operation (not needed when using copyobj)
end

xa = logspace(1,0,n);
xa = (xa-1); xa=xa/max(xa);

v = zeros(iData, n, 1);
for index=1:n
  this = a.*xa(index) + b.*(1-xa(index));
  this = setalias(this,'ratio',[ xa(index) 1-xa(index) ]);
  v(index) = this;
end

