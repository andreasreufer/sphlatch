function _out = bf_nbody( _in, _G, _dt)
  _out = _in;

  noDims  = size(_in.pos)(1);
  noParts = size(_in.pos)(2);

  # drift
  _out.oacc = _out.acc;
  _out.pos += _out.vel*_dt + 0.5*_dt*_dt*_out.acc;

  # brute-force nbody calculation
  _out.acc = zeros(noDims, noParts);
  for i = 1:noParts
    for j = 1:noParts
      if i != j
        r = _in.pos(:,i) - _in.pos(:,j);
        rdist = sqrt( sum( r .* r ) );

        _out.acc(:,i) -= (_G*_in.m(j)/(rdist^3.))*r;
      end
    end
  end

  # kick
  _out.vel += 0.5*(_out.acc + _out.oacc)*_dt;

endfunction


