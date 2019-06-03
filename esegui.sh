cd progetto/nasm32
nasm -f elf32 compute_residual_opt.nasm
cd ..
gcc -O0 -m32 -msse compute_residual_opt.o pqnn32c.c -o pqnn32c
./pqnn32c prova/prova -noexaustive -kc 100 -k 40 -nr 7000