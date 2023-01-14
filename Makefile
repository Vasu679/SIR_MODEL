SiRmodel: SiRmodel.o
	c++ -o SiRmodel SiRmodel.o -ltrapfpe -lpgplot -lcpgplot -lX11

SiRmodel.o: SiRmodel.cpp
	c++ SiRmodel.cpp -c
