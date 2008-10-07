#
# h5part_slice.m, Andreas Reufer, 06.07.2008
#
# _out = h5part_slice( _in, _min, _max )
#

function _out = h5part_slice( _in, _min, _max )
  
  keep = (
    _min(1) < _in.pos(1,:) &
    _min(2) < _in.pos(2,:) &
    _min(3) < _in.pos(3,:) &
    _max(1) > _in.pos(1,:) &
    _max(2) > _in.pos(2,:) &
    _max(3) > _in.pos(3,:)
    );

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

