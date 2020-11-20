# some macros
CC		= C:/Rtools/mingw64/bin/gcc.exe
CFLAGS		= -std=gnu99 -O2
R_C_HEADER	= C:/"Program Files"/R/R-4.0.2/include
DLL_BLAS	= C:/"Program Files"/R/R-4.0.2/bin/x64/Rblas.dll 
DLL_LAPACK	= C:/"Program Files"/R/R-4.0.2/bin/x64/Rlapack.dll 
DLL_R		= C:/"Program Files"/R/R-4.0.2/bin/x64/R.dll

# link (note PowerShell and cmd fail when lf + cr is used by '/') 
cbacon: cbacon.o wquantile.o 
	$(CC) -shared -o $@.dll $@.o $(DLL_BLAS) $(DLL_LAPACK) $(DLL_R) wquantile.o

# compile
cbacon.o: cbacon.c 
	$(CC) -Wall -pedantic -I $(R_C_HEADER) -c $*.c $(CFLAGS) 

# compile
wquantile.o: wquantile.c 
	$(CC) -Wall -pedantic -I $(R_C_HEADER) -c $*.c $(CFLAGS) 

# housekeeping (using make clean)
.PHONY: clean
clean: 
	rm cbacon.o wquantile.o 
