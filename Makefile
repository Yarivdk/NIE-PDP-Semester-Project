output: src/main.o src/MWCSolver.o src/dataFunctions.o
	g++ src/main.o src/MWCSolver.o src/dataFunctions.o -o output

src/main.o: src/main.cpp headers/dataFunctions.h headers/MWCSolver.h
	g++ -c src/main.cpp -o src/main.o -Iheaders

src/MWCSolver.o: src/MWCSolver.cpp headers/MWCSolver.h
	g++ -c src/MWCSolver.cpp -o src/MWCSolver.o -Iheaders

src/dataFunctions.o: src/dataFunctions.cpp headers/dataFunctions.h
	g++ -c src/dataFunctions.cpp -o src/dataFunctions.o -Iheaders

clean: 
	rm -f src/*.o output
