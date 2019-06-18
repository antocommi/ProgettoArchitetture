cd progetto/nasm32
nasm -f elf32 distanza.nasm
nasm -f elf32 calcolaIndice.nasm
nasm -f elf32 calcolaPQ.nasm
nasm -f elf32 dist.nasm

cd ..
# gcc -O0 -m32 -msse  ./nasm32/calcolaIndice.o ./nasm32/calcolaPQ.o ./nasm32/compute_residual_opt.o ./nasm32/distanza.o pqnn32c_backup.c -o pqnn32c -lm
# ./"pqnn32c" prova/prova -noexaustive -kc 350 -nr 400 -adc -knn 1

# VALGRIND
gcc -g -O1 -m64 -no-pie -mavx ./nasm64/calcolaIndice.o ./nasm64/distanza.o ./nasm64/dist_matrix.o ./nasm64/dist.o ./nasm64/calcolaFob.o ./nasm64/calcolaPQ.o ./nasm64/calcolaSimmetriche.o ./nasm64/somma.o ./nasm64/creaMatriceDistanze.o ./pqnn64c.c -o ./pqnn64c -lm
man valgrind
valgrind --leak-check=yes -v ./pqnn64c ./prova/prova -noexaustive -kc 350 -nr 400 -knn 4
cd ..
#gcc -O3 -m32 ./progetto/pqnn32c.c -o ./progetto/pqnn32c -lm
#./progetto/pqnn32c ./progetto/prova/prova -noexaustive -kc 350 -nr 400 -adc
