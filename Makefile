sample.out:	main.o init.o insert.o update.o parameter.o print.o 
		gcc -o sample.out ./object/main.o ./object/init.o ./object/insert.o ./object/update.o ./object/print.o ./object/parameter.o -lm -fopenmp -O3
main.o:	./main.c
		gcc -o ./object/main.o -c ./main.c
init.o: ./src/init.c
		gcc -o ./object/init.o -c ./src/init.c
insert.o: ./src/insert.c
		gcc -o ./object/insert.o -c ./src/insert.c -lm
update.o: ./src/update.c
		gcc -o ./object/update.o -c ./src/update.c -lm -fopenmp -O3
print.o: ./src/print.c
		gcc -o ./object/print.o -c ./src/print.c
parameter.o: ./src/parameter.c
		gcc -o ./object/parameter.o -c ./src/parameter.c