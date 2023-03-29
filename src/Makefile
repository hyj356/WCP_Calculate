FC = gfortran
FILES = mod_type.f90 wcsro.f90
PARRALLEL = -O2 -fopenmp -fmax-errors=1 -Wall -Wno-unused-function -Wno-maybe-uninitialized -Wno-unused-variable
SERIAL = -O2 -fmax-errors=1 -Wall -Wno-unused-function -Wno-maybe-uninitialized -Wno-unused-variable

.PHONY:clean

wcp_op: ${FILES}	# 采用openMP并行编译
	${FC} ${FILES} -o $@ ${PARRALLEL}

wcp_se: ${FILES}	# 采用串行程序
	${FC} ${FILES} -o $@ ${SERIAL}

clean:
	rm -rf *.mod