main: main.o mylib.o support.o MNA.o delay.o STA.o sizing.o EPA.o delaycalc.o
	g++ -o asloptimizer main.o mylib.o support.o MNA.o delay.o STA.o sizing.o EPA.o delaycalc.o -std=c++11

delaycalc.o: delaycalc.cpp gate.h
	g++ -c delaycalc.cpp -std=c++11
EPA.o: EPA.cpp gate.h
	g++ -c EPA.cpp -std=c++11
sizing.o: sizing.cpp gate.h
	g++ -c sizing.cpp -std=c++11
STA.o: STA.cpp gate.h
	g++ -c STA.cpp -std=c++11
delay.o: delay.cpp gate.h
	g++ -c delay.cpp -std=c++11
MNA.o: MNA.cpp gate.h
	g++ -c MNA.cpp -std=c++11
support.o: support.cpp gate.h
	g++ -c support.cpp -std=c++11
mylib.o: mylib.cpp gate.h
	g++ -c mylib.cpp -std=c++11
main.o: main.cpp gate.h
	g++ -c main.cpp -std=c++11

clean:
	rm -f main main.o mylib.o support.o MNA.o delay.o STA.o sizing.o EPA.o delaycalc.o
