#
# h5part_slice.m, Andreas Reufer, 23.07.2008
#
# _idx = h5part_max( _dump, _key, _n )
#

function _idx = h5part_max( _dump, _key, _n )
  [x, idx] = sort( getfield( _dump, _key ), 'descend' );
  _idx = idx(1:_n);
endfunction

