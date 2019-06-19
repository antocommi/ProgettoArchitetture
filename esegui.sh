# cd progetto/nasm32
# nasm -f elf32 distanza.nasm
# nasm -f elf32 calcolaIndice.nasm
# nasm -f elf32 calcolaPQ.nasm
# nasm -f elf32 dist.nasm

# cd ..
# gcc -O0 -m32 -msse  ./nasm32/calcolaIndice.o ./nasm32/calcolaPQ.o ./nasm32/compute_residual_opt.o ./nasm32/distanza.o pqnn32c_backup.c -o pqnn32c -lm
# ./"pqnn32c" prova/prova -noexaustive -kc 350 -nr 400 -adc -knn 1
# cd ..
# VALGRIND
cd progetto/nasm64
nasm -f elf64 calcolaIndice.nasm
nasm -f elf64 distanza.nasm
nasm -f elf64 dist_matrix.nasm
nasm -f elf64 dist.nasm
nasm -f elf64 calcolaFob.nasm
nasm -f elf64 calcolaPQ.nasm
nasm -f elf64 calcolaSimmetriche.nasm
nasm -f elf64 somma.nasm
nasm -f elf64 creaMatriceDistanze.nasm
nasm -f elf64 compute_residual.nasm
cd ..
gcc -O1 -m64 -g -no-pie -mavx ./nasm64/calcolaIndice.o ./nasm64/distanza.o ./nasm64/dist_matrix.o ./nasm64/dist.o ./nasm64/calcolaFob.o ./nasm64/calcolaPQ.o ./nasm64/calcolaSimmetriche.o ./nasm64/somma.o ./nasm64/creaMatriceDistanze.o ./nasm64/compute_residual.o pqnn64c.c -o pqnn64c
valgrind --leak-check=yes --read-inline-info=yes --track-origins=yes --read-var-info=yes -v ./pqnn64c ./prova/prova -noexaustive -kc 350 -nr 400 -knn 1 -sdc
cd ..
# gcc -g -O1 -m64 -no-pie -mavx ./nasm64/calcolaIndice.o ./nasm64/distanza.o ./nasm64/dist_matrix.o ./nasm64/dist.o ./nasm64/calcolaFob.o ./nasm64/calcolaPQ.o ./nasm64/calcolaSimmetriche.o ./nasm64/somma.o ./nasm64/creaMatriceDistanze.o ./pqnn64c.c -o ./pqnn64c -lm
# man valgrind
# valgrind --leak-check=yes -v ./pqnn64c ./prova/prova -noexaustive -kc 350 -nr 400 -knn 4
# cd ..
#gcc -O3 -m32 ./progetto/pqnn32c.c -o ./progetto/pqnn32c -lm
#./progetto/pqnn32c ./progetto/prova/prova -noexaustive -kc 350 -nr 400 -adc
