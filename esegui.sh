cd progetto/nasm32
nasm -f elf32 distanza.nasm
cd ..
gcc -O0 -m32 -msse ./nasm32/distanza.o pqnn32c.c -o pqnn32c -lm
./"pqnn32c" prova/prova -noexaustive -kc 100 -k 40 -nr 7000
cd ..

# VALGRIND
#gcc -g -O1 -m32 ./progetto/pqnn32c.c -o ./progetto/pqnn32c
#valgrind --leak-check=yes -v ./progetto/pqnn32c ./progetto/prova/prova -noexaustive -kc 350 -nr 400 -knn 4

#gcc -O3 -m32 ./progetto/pqnn32c.c -o ./progetto/pqnn32c -lm
#./progetto/pqnn32c ./progetto/prova/prova -noexaustive -kc 350 -nr 400 -adc
