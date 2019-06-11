# cd progetto/nasm32
# nasm -f elf32 compute_residual_opt.nasm
# cd ..
# gcc -O0 -lm -m32 -msse ./nasm32/compute_residual_opt.o pqnn32c.c -o pqnn32c
# ./"pqnn32c" prova/prova -noexaustive -kc 100 -k 40 -nr 7000
# cd ..

# VALGRIND
gcc -g -O1 -m32 ./progetto/pqnn32c.c -o ./progetto/pqnn32c2
valgrind --leak-check=yes -v ./progetto/pqnn32c2 ./progetto/prova/prova -noexaustive -kc 100 -nr 400
