F2PY3 = f2py3
F2PY2 = f2py

PY2 := $(shell command -v $(F2PY2) 2> /dev/null)
PY3 := $(shell command -v $(F2PY3) 2> /dev/null)

all:
ifdef PY2
	$(MAKE) ossir.so
endif
ifdef PY3
	$(MAKE) ossir.cpython-34m.so
endif

ossir.cpython-34m.so: oss_ir.F90
	$(F2PY3) -c -m --f90flags="-ffast-math -Wno-tabs -ffree-form" ossir $<

ossir.so: oss_ir.F90
	$(F2PY2) -c -m --f90flags="-ffast-math -Wno-tabs -ffree-form" ossir $<

clean:
	rm -f ossir*.so
