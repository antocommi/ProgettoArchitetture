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

; %include "sseutils.nasm"

section .data			; Sezione contenente dati inizializzati

uno:		dd		1.0
;
;align 16
;inizio:		dd		1.0, 2.0, 3.0, 4.0

section .bss			; Sezione contenente dati non inizializzati

;alignb 16
;vec2:		resq	4

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

global assegnaSoluzioneF

pag		equ		8
xk              equ             12
n               equ             16

assegnaSoluzioneF:
		push		ebp; salva il Base Pointer
		mov		ebp, esp; il Base Pointer punta al Record di Attivazione corrente
		push		ebx; salva i registri da preservare
		push		esi
		push		edi
		
		mov		eax, [ebp+pag]
                mov             ebx, [ebp+xk]
                mov             ecx, [ebp+n]

                xor esi, esi
          ciclo:movaps xmm0, [ebx+esi*4]
                movaps xmm1, [ebx+esi*4+16]
                movaps xmm2, [ebx+esi*4+32]
                movaps xmm3, [ebx+esi*4+48]
                movaps xmm4, [ebx+esi*4+64]
                movaps xmm5, [ebx+esi*4+80]
                movaps xmm6, [ebx+esi*4+96]
                movaps xmm7, [ebx+esi*4+112]
                movaps [eax+esi*4], xmm0
                movaps [eax+esi*4+16], xmm1
                movaps [eax+esi*4+32], xmm2
                movaps [eax+esi*4+48], xmm3
                movaps [eax+esi*4+64], xmm4
                movaps [eax+esi*4+80], xmm5
                movaps [eax+esi*4+96], xmm6
                movaps [eax+esi*4+112], xmm7
                xorps xmm0, xmm0
                movaps [ebx+esi*4], xmm0
                movaps [ebx+esi*4+16], xmm0
                movaps [ebx+esi*4+32], xmm0
                movaps [ebx+esi*4+48], xmm0
                movaps [ebx+esi*4+64], xmm0
                movaps [ebx+esi*4+80], xmm0
                movaps [ebx+esi*4+96], xmm0
                movaps [ebx+esi*4+112], xmm0
                add esi, 32
                cmp esi, ecx
                jl ciclo

		pop	edi; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp; ripristina lo Stack Pointer
		pop	ebp; ripristina il Base Pointer
		ret	   ; torna alla funzione C chiamante
