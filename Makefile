all: ising cyclicPerturbation cyclicPerturbationFast cyclicPerturbationKrylov cyclicPerturbationKrylov2 cyclicPerturbationTest test spectralFunctionKrylov spectralFunctionKrylov2

ising: src/ising.cpp
	g++ -c ./src/ising.cpp -fopenmp -o ./lib/ising.o -O3 -I./include -I./../Libraries/eigen

cyclicPerturbation: src/cyclicPerturbation.cpp
	g++ -c ./src/cyclicPerturbation.cpp -fopenmp -o ./lib/cyclicPerturbation.o -O3 -I./include -I./../Libraries/eigen -I./../Libraries/json/include
	g++ -o _cyclicPerturbation.out -fopenmp -O3 ./lib/cyclicPerturbation.o ./lib/ising.o

cyclicPerturbationFast: src/cyclicPerturbationFast.cpp
	g++ -c ./src/cyclicPerturbationFast.cpp -fopenmp -o ./lib/cyclicPerturbationFast.o -O3 -I./include -I./../Libraries/eigen -I./../Libraries/json/include
	g++ -o _cyclicPerturbationFast.out -fopenmp -O3 ./lib/cyclicPerturbationFast.o ./lib/ising.o

cyclicPerturbationKrylov: src/cyclicPerturbationKrylov.cpp
	g++ -c ./src/cyclicPerturbationKrylov.cpp -fopenmp -o ./lib/cyclicPerturbationKrylov.o -O3 -I./include -I./../Libraries/eigen -I./../Libraries/json/include
	g++ -o _cyclicPerturbationKrylov.out -fopenmp -O3 ./lib/cyclicPerturbationKrylov.o ./lib/ising.o

cyclicPerturbationKrylov2: src/cyclicPerturbationKrylov2.cpp
	g++ -c ./src/cyclicPerturbationKrylov2.cpp -fopenmp -o ./lib/cyclicPerturbationKrylov2.o -O3 -I./include -I./../Libraries/eigen -I./../Libraries/json/include
	g++ -o _cyclicPerturbationKrylov2.out -fopenmp -O3 ./lib/cyclicPerturbationKrylov2.o ./lib/ising.o

spectralFunctionKrylov: src/spectralFunctionKrylov.cpp
	g++ -c ./src/spectralFunctionKrylov.cpp -fopenmp -o ./lib/spectralFunctionKrylov.o -O3 -I./include -I./../Libraries/eigen -I./../Libraries/json/include
	g++ -o _spectralFunctionKrylov.out -fopenmp -O3 ./lib/spectralFunctionKrylov.o ./lib/ising.o

spectralFunctionKrylov2: src/spectralFunctionKrylov2.cpp
	g++ -c ./src/spectralFunctionKrylov2.cpp -fopenmp -o ./lib/spectralFunctionKrylov2.o -O3 -I./include -I./../Libraries/eigen -I./../Libraries/json/include
	g++ -o _spectralFunctionKrylov2.out -fopenmp -O3 ./lib/spectralFunctionKrylov2.o ./lib/ising.o

cyclicPerturbationTest: src/cyclicPerturbationTest.cpp
	g++ -c ./src/cyclicPerturbationTest.cpp -fopenmp -o ./lib/cyclicPerturbationTest.o -O3 -I./include -I./../Libraries/eigen -I./../Libraries/json/include
	g++ -o _cyclicPerturbationTest.out -fopenmp -O3 ./lib/cyclicPerturbationTest.o ./lib/ising.o

test: src/test.cpp
	g++ -c ./src/test.cpp -fopenmp -o ./lib/test.o -O3 -I./include -I./../Libraries/eigen -I./../Libraries/json/include
	g++ -o _test.out -fopenmp -O3 ./lib/test.o ./lib/ising.o