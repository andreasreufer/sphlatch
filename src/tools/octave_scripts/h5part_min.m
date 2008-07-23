#
# h5part_slice.m, Andreas Reufer, 23.07.2008
#
# _idx = h5part_min( _dump, _key, _n )
#

function _idx = h5part_min( _dump, _key, _n )

  [x, idx] = sort( getfield( _dump, _key ), 'ascend' );

  printf("  ID        %s\n", _key );
  printf("  %i        %e\n", _dump.id( idx(1:_n) ), x(1:_n ) );

  _idx = idx(1:_n);
endfunction

