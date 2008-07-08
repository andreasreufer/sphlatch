#
# h5part_compare.m, Andreas Reufer, 06.07.2008
#
# delta = h5part_compare( ref, cmp )
#
# compare a reference particle structure ("ref") to
# another particle structure ("cmp") and store the
# difference delta = ref - cmp
#
# it assumes the same set of particles in both structures,
# namely the same set of particle ids!
#

function delta = h5part_compare(ref, cmp)

  # sort both ref and cmp ids
  [ ref.id_sorted, ref.idx_sorted] = sort( ref.id );
  [ cmp.id_sorted, cmp.idx_sorted] = sort( cmp.id );

  # create empty delta structure
  delta = struct;

  delta.idx_ref = ref.idx_sorted;
  delta.idx_cmp = cmp.idx_sorted;

  for [refVal, refKey] = ref

    # check whether both ref and cmp have the variable with the current key
    useKey = true;
    try
      getfield( cmp, refKey );
    catch
      useKey = false;
    end

    if ( strcmp( refKey, "id_sorted") || strcmp( refKey, "idx_sorted") )
      useKey = false;
    end

    if ( useKey )
      printf("delta of %s\n", refKey);
      delta = setfield( delta, strcat("d_", refKey),
        refVal(:,ref.idx_sorted) -
        (getfield( cmp, refKey ))(:,cmp.idx_sorted) );
    endif

  endfor

  # also store the reference particle positions
  delta = setfield( delta, "refpos", ref.pos(:, ref.idx_sorted) );

endfunction

