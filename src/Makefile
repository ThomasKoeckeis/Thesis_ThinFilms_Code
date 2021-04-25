# MAKEFILE

CC	= gcc
FLAGS	= -g -Wall -pedantic
INCLS	= -I/.
LIBS	= -lm

.c.o:  ; $(CC) -c $(FLAGS) $(INCL) $*.c

OBJ_NUC=	main.o \
		allocate.o \
		nucleate.o \
		intersection.o \
		boundary.o \
		contact.o \
		stress.o \
		graphics.o \
		auxiliary.o

grow:	$(OBJ_NUC)
	$(CC) $(FLAGS) $(INCL) $(OBJ_NUC) -o grow $(LIBS)

clean:
	rm -f *.o core
