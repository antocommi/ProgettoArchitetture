; ---------------------------------------------------------
; PageRank con istruzioni SSE a 32 bit
; ---------------------------------------------------------
; F. Angiulli
; 23/11/2017
;

;
; Software necessario per l'esecuzione:
;
;     NASM (www.nasm.us)
;     GCC (gcc.gnu.org)
;
; entrambi sono disponibili come pacchetti software 
; installabili mediante il packaging tool del sistema 
; operativo; per esempio, su Ubuntu, mediante i comandi:
;
;     sudo apt-get install nasm
;     sudo apt-get install gcc
;
; potrebbe essere necessario installare le seguenti librerie:
;
;     sudo apt-get install lib32gcc-4.8-dev (o altra versione)
;     sudo apt-get install libc6-dev-i386
;
; Per generare file oggetto:
;
;     nasm -f elf32 abod32.nasm 
;

 %include "sseutils.nasm"

section .data			; Sezione contenente dati inizializzati
;align 16
uno:		dq		1.0
;
;align 16
;inizio:		dd		1.0, 2.0, 3.0, 4.0

section .bss			; Sezione contenente dati non inizializzati

;alignb 16
somma:		resq	1
ri:             resd    1
ri2:            resq    1

section .text			; Sezione contenente il codice macchina


; ----------------------------------------------------------
; macro per l'allocazione dinamica della memoria
;
;	getmem	<size>,<elements>
;
; alloca un'area di memoria di <size>*<elements> bytes
; (allineata a 16 bytes) e restituisce in EAX
; l'indirizzo del primo bytes del blocco allocato
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)
;
;	fremem	<address>
;
; dealloca l'area di memoria che ha inizio dall'indirizzo
; <address> precedentemente allocata con getmem
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)

extern get_block
extern free_block
extern sommad

%macro	getmem	2
	mov	eax, %1
	push	eax
	mov	eax, %2
	push	eax
	call	get_block
	add	esp, 8
%endmacro

%macro	fremem	1
	push	%1
	call	free_block
	add	esp, 4
%endmacro


; ------------------------------------------------------------
; Funzione pagerank32
; ------------------------------------------------------------

global prodSparseF

row		equ		8
col             equ             12
vd              equ             16
N               equ             20
xk1             equ             24
xk              equ             28
d               equ             32
M               equ             36

prodSparseF:
		push		ebp; salva il Base Pointer
		mov		ebp, esp; il Base Pointer punta al Record di Attivazione corrente
		push		ebx; salva i registri da preservare
		push		esi
		push		edi

                xor esi, esi
                mov edx,[ebp+xk1]
                mov ebx,[ebp+d]
                mov eax,[ebp+vd]
                mov edi,[ebp+M]
         cicloV:movaps xmm0,[edx+esi*4]
                movaps xmm1,[edx+esi*4+16]
                movaps xmm2,[edx+esi*4+32]
                movaps xmm3,[edx+esi*4+48]
                movaps xmm4,[edx+esi*4+64]
                movaps xmm5,[edx+esi*4+80]
                movaps xmm6,[edx+esi*4+96]
                movaps xmm7,[edx+esi*4+112]
                mulps xmm0,[ebx+esi*4]
                mulps xmm1,[ebx+esi*4+16]
                mulps xmm2,[ebx+esi*4+32]
                mulps xmm3,[ebx+esi*4+48]
                mulps xmm4,[ebx+esi*4+64]
                mulps xmm5,[ebx+esi*4+80]
                mulps xmm6,[ebx+esi*4+96]
                mulps xmm7,[ebx+esi*4+112]
                movaps [eax+esi*4],xmm0
                movaps [eax+esi*4+16],xmm1
                movaps [eax+esi*4+32],xmm2
                movaps [eax+esi*4+48],xmm3
                movaps [eax+esi*4+64],xmm4
                movaps [eax+esi*4+80],xmm5
                movaps [eax+esi*4+96],xmm6
                movaps [eax+esi*4+112],xmm7
                add esi, 32
                cmp esi,edi
                jl cicloV
                xor esi,esi
                xor edi, edi
                mov edx, [ebp+row]
                mov ebx,[ebp+vd]
                mov ecx, [ebp+col]
                mov eax, [ecx]
          ciclo:
                xorps xmm2,xmm2
                cmp eax, 0
                je fine
                
         ciclor:
                mov ecx,[edx+edi*4]
                sal ecx,2 
                add ecx,ebx
                movss xmm0,[ecx]
                addss xmm2,xmm0
   
                inc edi
                dec eax
                cmp eax, 0
                jg ciclor
                
           fine:mov ecx,[ebp+xk]
                movss [ecx+esi*4],xmm2
                mov ecx, [ebp+col]
                mov eax, [ecx+esi*4+4]
                sub eax, [ecx+esi*4]
                inc esi  
                mov ecx, [ebp+N]
                cmp esi, ecx
                jl ciclo
                

		pop	edi; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp; ripristina lo Stack Pointer
		pop	ebp; ripristina il Base Pointer
		ret	   ; torna alla funzione C chiamante
