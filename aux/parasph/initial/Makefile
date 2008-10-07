all:
	$(CC) $(CFLAGS) -c xdrlib.c
	$(FC) $(FFLAGS) -c libaneos.f
	$(FC) $(FFLAGS) -lm -lgfortran -lm \
	  xdrlib.o libaneos.o setup-collision.f \
	  -o setup-collision

clean:
	rm libaneos.o xdrlib.o setup-collision
