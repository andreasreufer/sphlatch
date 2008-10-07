#
# h5part_projvect.m, Andreas Reufer, 14.08.2008
#
# _dist = h5part_projvect( _in, _vect, _center)
#

function _out = h5part_projvect( _in, _vect, _center = [0.,0.,0.])
  relpos(1,:) = _in.pos(1,:) - _center(1);
  relpos(2,:) = _in.pos(2,:) - _center(2);
  relpos(3,:) = _in.pos(3,:) - _center(3);
  dist = sqrt( relpos(1,:).*relpos(1,:) +
               relpos(2,:).*relpos(2,:) +
               relpos(3,:).*relpos(3,:) );

  _in = setfield( _in, strcat(_vect, "_proj"),
            sum( (getfield(_in, _vect)(:,:)) .* (_in.pos(:,:)) ) ./ dist );
  
  _in = setfield( _in, "dist", dist );
  
  _out = _in;
endfunction

