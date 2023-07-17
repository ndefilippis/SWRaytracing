function fa = avg(f,d,periodic,shift,endoff)

%  f = avg(f,dim,periodic,shift,endoff) Average 2 or 3 dim field f along dimension
%  d (=1 or 2).  Periodic flag is optionl, with default false.  If false,
%  fa(end) = f(end)/2.  If true, fa(end) = (f(end)+f(1))/2. Optional
%  flag shift (default false) shifts result so that periodic element
%  is in first position of fa instead of last. Optional endoff
%  (default false) truncates the last element along the averaging dimension.

if (nargin<3), endoff=false; shift=false; periodic=false; end
if (nargin<4), endoff=false; shift=false; end
if (nargin<5), endoff=false; end
fa = zeros(size(f));

switch d
  case 1
    fa(1:end-1,:,:) = (f(1:end-1,:,:) + f(2:end,:,:))/2;
    fa(end,:,:) = f(end,:,:)/2;
    if (periodic), fa(end,:,:) = (f(1,:,:)+f(end,:,:))/2; end 
    if (shift),    fa = circshift(fa,[1 0 0]); end
    if (endoff),   fa = fa(1:end-1,:,:);  end
  case 2
    fa(:,1:end-1,:) = (f(:,1:end-1,:) + f(:,2:end,:))/2;
    fa(:,end,:) = f(:,end,:)/2;
    if (periodic), fa(:,end,:) = (f(:,1,:)+f(:,end,:))/2; end 
    if (shift),    fa = circshift(fa,[0 1 0]); end
    if (endoff),   fa = fa(:,1:end-1,:);  end
  otherwise
    message('Averaging dimension unrecognized in avg:',d)
end

