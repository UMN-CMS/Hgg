photosLHEF: photosLHEF.o photos.o
	g++ -o photosLHEF $^ -lgfortran

photosLHEF.o: photosLHEF.cc
	gcc -c $^ -Wall -O2 

photos.o: photos.F
	gfortran -c photos.F -fno-automatic -g -ffixed-line-length-132 -fno-backslash  -DCZTERYD
