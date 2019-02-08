CC = gcc
LFLAGS = -lflint -lgmp -lm

%.o: %.c
	$(CC) -fPIC -c $<

libembed.so: AE.o h90.o minpoly.o tensor.o basis_change.o nth-root.o linfactor.o import_conway.o
	gcc -shared $+ $(LFLAGS) -o libembed.so

AE.test.exe: AE.o test_AE.o

h90.test.exe: h90.o AE.o test_h90.o

minpoly.test.exe: minpoly.o test_minpoly.o

tensor.test.exe: tensor.o test_tensor.o

tensor_scalarmul.test.exe: tensor.o test_tensor_scalarmul.o

basis_change.test.exe: basis_change.o minpoly.o test_basis_change.o

%.test.exe:
	$(CC) $+ -o $@ $(LFLAGS)

%.test: %.test.exe
	./$<

all-tests := $(addsuffix .test,$(subst test_,,$(basename $(wildcard test_*.c))))
test: $(all-tests)

.PHONY: test
