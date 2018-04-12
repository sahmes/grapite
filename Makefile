GRAPELIB=../yebisu/libyebisug6.a
ETICSLIB=../etics/src/libetics.a
NVCC=nvcc -ccbin g++-6
GPU_ARCH = 61

GRAPELIBFULLPATH=$(shell readlink -f $(GRAPELIB))
ETICSLIBFULLPATH=$(shell readlink -f $(ETICSLIB))
default: libgrapite.a

libgrapite.a:
	$(NVCC) -arch=sm_$(GPU_ARCH) -O3 -c grapite.cu
	rm -rf grapelibtmp
	mkdir grapelibtmp
	cd grapelibtmp; ar x $(GRAPELIBFULLPATH); \
	for i in $$(ls); do \
		objcopy --redefine-syms=../redefine-symbols.txt $$i prefixed_$$i; \
		rm -f $$i; \
	done; ar x $(ETICSLIBFULLPATH)
	ar -r libgrapite.a grapite.o grapelibtmp/*.o
	rm -rf grapelibtmp
	ranlib libgrapite.a

	
test: libgrapite.a
	cc -o test test.c -fopenmp -L. -L/usr/local/cuda/lib64 -lgrapite -lcuda -lcudart
clean:
	rm -f *.o *.a test
	rm -rf grapelibtmp
