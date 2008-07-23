#
# h5part_find.m, Andreas Reufer, 23.07.2008
#
# _idx = h5part_findID( _dump, _id )
#

function _idx = h5part_findID( _dump, _id)
  [ids, _idx] = find( round( _dump.id ) == round( _id ) );
endfunction

