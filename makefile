all: sphere

sphere: main.o kVector.o spherelit.o animation.o
	g++ main.o kVector.o spherelit.o animation.o -o sphere -lSDL -lm `sdl-config --cflags --libs`

main.o: main.cc
	g++ -c main.cc

animation.o: animation.cc
	g++ -c animation.cc

spherelit.o: spherelit.cc
	g++ -c spherelit.cc

kVector.o: kVector.cc
	g++ -c kVector.cc
