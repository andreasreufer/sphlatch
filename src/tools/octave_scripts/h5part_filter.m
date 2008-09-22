#
# h5part_slice.m, Andreas Reufer, 06.07.2008
#
# _out = h5part_filter( _in, _smaller, _bigger )
#

function _out = h5part_filter( _in, _smaller, _bigger )
  
  keep = ( _smaller < _bigger );

  # build a vector of indices to be kept
  noParts = size( _in.pos )(2);
  curKeep = 1;
  for i = 1:noParts
    if keep(i)
      keepIdx( curKeep ) = i;
      curKeep++;
    end
  end

  printf("keeping %i of %i particles\n", sum(keep), noParts );

  # create the _out structure
  _out = struct;
  for [_inVal, _inKey] = _in
    _out = setfield( _out, _inKey, (getfield( _in, _inKey))(:,keepIdx) );
  endfor

endfunction

