# turn object file into merge executable (linking stage)
# -g for debugging info
nbod: main.o
	g++ -g -o nbod main.o

# create object file from source (compilation stage)
main.o: main.cpp
	g++ -c main.cpp

# remove object file and executable
clean:
	rm *.o nbod