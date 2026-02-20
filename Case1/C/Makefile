CC = gcc
CFLAGS = -O2 -I/usr/include/lapacke
LIBS = -llapacke -llapack -lblas -lm

all:
	$(CC) $(CFLAGS) Case1.c TSQR.c -o Case1 $(LIBS)
	$(CC) $(CFLAGS) Scaling.c TSQR.c -o Scaling $(LIBS)

clean:
	rm -f Case1 TSQR Scaling
