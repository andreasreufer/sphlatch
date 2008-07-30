#
# h5part_slice.m, Andreas Reufer, 06.07.2008
#
# _dist = h5part_dist( _in)
#

function _dist = h5part_dist( _in, _center )

  _dist = sqrt(
            ( _in.pos(1,:) - _center(1) ).*( _in.pos(1,:) - _center(1) ) +
            ( _in.pos(2,:) - _center(2) ).*( _in.pos(2,:) - _center(2) ) +
            ( _in.pos(3,:) - _center(3) ).*( _in.pos(3,:) - _center(3) ) );

endfunction

