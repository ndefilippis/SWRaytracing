function movie_out = make_swmovie(F,field,ptype,cmap,xvec,yvec,cvec);

% [movie_out] = MAKE_SWMOVIE(field,ptype,cmap,xvec,yvec,cvec)  
%     Make movie of 2d 'field' w/ PLOTFIELD.  Other args optional:
%     OPTION      FUNCTION             DEFAULT    
%     'ptype'     Plot type            'pc' (pcolor w/ colorbar)
%     'cmap'      Colormap             'jet'
%     xvec        x coordinate vector  1:size(field,1)
%     yvec        y coordinate vector  1:size(field,2)
%     cvec        color axis value     auto
%
%  This version assumes F is a structure, and field is the fieldname.
%
%     See also PLOTFIELD, MOVIEIN, GETFRAME, MAKE_FILEMOVIE.

% Defaults:
cvecd = 'auto'; yvecd = 1:size(field,2); xvecd = 1:size(field,1);
cmapd = 'default'; ptyped = 'pc';

switch nargin
  case 1, cvec = cvecd; yvec = yvecd; xvec = xvecd; cmap = cmapd;
          ptype = ptyped;    
  case 2, cvec = cvecd; yvec = yvecd; xvec = xvecd; cmap = cmapd; 
  case 3, cvec = cvecd; yvec = yvecd; xvec = xvecd;  
  case 4, cvec = cvecd; yvec = yvecd;
  case 5, cvec = cvecd;
end

movie_out = moviein(size(F));  % pre-allocate memory
clf

for frm=1:size(field,3)
  disp(frm)
  f=F.field(:,:,frm);
  plotfield(f,0,ptype,cmap,xvec,yvec,cvec); 
  movie_out(:,frm) = getframe;
end
