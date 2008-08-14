#
# h5part_slice.m, Andreas Reufer, 14.08.2008
#
# _dist = h5part_projvect( _in, _vect, _center)
#

function _out = h5part_projvect( _in, _vect, _center = [0.,0.,0.])
  dist = sqrt(
            ( _in.pos(1,:) - _center(1) ).*( _in.pos(1,:) - _center(1) ) +
            ( _in.pos(2,:) - _center(2) ).*( _in.pos(2,:) - _center(2) ) +
            ( _in.pos(3,:) - _center(3) ).*( _in.pos(3,:) - _center(3) ) );
  _out = _in;
endfunction

