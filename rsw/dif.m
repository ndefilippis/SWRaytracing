function fd = dif(f,d,periodic,shift,endoff)

%  fd = dif(f,dim,periodic,shift) Difference 2 dim field f along
%  dimension d.  Periodic flag is optional, with default false.  If
%  false, fd(end) = -f(end).  If true, fd(end) = f(1)-f(end). Optional flag
%  shift (default false) shifts result so that periodic element is in
%  first position of fd instead of last.

if (nargin<3), endoff=false; shift=false; periodic=false; end
if (nargin<4), endoff=false; shift=false; end
if (nargin<5), endoff=false; end
fd = zeros(size(f));

switch d
  case 1
    fd(1:end-1,:,:) = f(2:end,:,:) - f(1:end-1,:,:);
    fd(end,:,:) = -f(end,:,:);
    if (periodic), fd(end,:,:) = f(1,:,:)-f(end,:,:); end
    if (shift),    fd = circshift(fd,[1 0 0]); end
    if (endoff),   fd = fd(1:end-1,:,:);  end
  case 2
    fd(:,1:end-1,:) = f(:,2:end,:) - f(:,1:end-1,:);
    fd(:,end,:) = -f(:,end,:);
    if (periodic), fd(:,end,:) = f(:,1,:)-f(:,end,:); end 
    if (shift),    fd = circshift(fd,[0 1 0]); end
    if (endoff),   fd = fd(:,1:end-1,:);  end
  otherwise
    message('Differencing dimension unrecognized in dif:',d)
end

