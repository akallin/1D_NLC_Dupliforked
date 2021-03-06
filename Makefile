CC = g++
OBJS = NLC_2D_TFIM.o CPU/magnetization.o CPU/GenHam.o CPU/Lanczos_07.o graphs.o
CFLAGS = -O2 -fopenmp -Wall -Wextra --pedantic

2d-cpu.out: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o 2d-cpu.out $(HEISLIBS) $(LANCLIBS)

NLC_2D_TFIM.o : NLC_2D_TFIM.cpp CPU/GenHam.h CPU/Lanczos_07.h CPU/simparam.h ../../BossBranch/Graphs/graphs.h
	$(CC) $(CFLAGS) -c NLC_2D_TFIM.cpp 

CPU/magnetization.o: CPU/magnetization.cpp CPU/magnetization.h
	$(CC) $(CFLAGS) -c CPU/magnetization.cpp -o CPU/magnetization.o

CPU/GenHam.o: CPU/GenHam.cpp CPU/GenHam.h CPU/Lanczos_07.h
	$(CC) $(CFLAGS) -c CPU/GenHam.cpp -o CPU/GenHam.o

CPU/Lanczos_07.o: CPU/Lanczos_07.cpp CPU/GenHam.h CPU/Lanczos_07.h
	$(CC) $(CFLAGS) -c CPU/Lanczos_07.cpp -o CPU/Lanczos_07.o

graphs.o: ../../BossBranch/Graphs/graphs.cpp ../../BossBranch/Graphs/graphs.h
	$(CC) $(CFLAGS) -c ../../BossBranch/Graphs/graphs.cpp

clean :
	rm CPU/*.o && rm *.o
