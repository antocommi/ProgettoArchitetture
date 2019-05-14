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

;uno:		dd		1.0
;
;align 16
;inizio:		dd		1.0, 2.0, 3.0, 4.0

section .bss			; Sezione contenente dati non inizializzati

;alignb 16
ri:             resd 1
ri2: 		resd	1

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

global sommaf

v		equ		8
dim             equ             12
s               equ             16

sommaf:
		push		ebp; salva il Base Pointer
		mov		ebp, esp; il Base Pointer punta al Record di Attivazione corrente
		push		ebx; salva i registri da preservare
		push		esi
		push		edi
		
		mov		eax, [ebp+v]
                mov             ebx, [ebp+dim]
                mov             ecx, [ebp+s]
 
                movaps xmm0, [eax]
                movaps xmm1, [eax+16]
                movaps xmm2, [eax+32]
                movaps xmm3, [eax+48]
                movaps xmm4, [eax+64]
                movaps xmm5, [eax+80]
                movaps xmm6, [eax+96]
                movaps xmm7, [eax+112]
                cmp ebx, 32
                je fine
                mov esi, 32
         cicloq:addps xmm0, [eax+esi*4]
                addps xmm1, [eax+esi*4+16]
                addps xmm2, [eax+esi*4+32]
                addps xmm3, [eax+esi*4+48]
                addps xmm4, [eax+esi*4+64]
                addps xmm5, [eax+esi*4+80]
                addps xmm6, [eax+esi*4+96]
                addps xmm7, [eax+esi*4+112]
                add esi, 32
                cmp esi,ebx
                jl cicloq 
           fine:addps xmm0, xmm1
                addps xmm0, xmm2
                addps xmm0, xmm3
                addps xmm0, xmm4
                addps xmm0, xmm5
                addps xmm0, xmm6
                addps xmm0, xmm7
                haddps xmm0, xmm0
                haddps xmm0, xmm0
                movss [ecx], xmm0	
                

		pop	edi; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp; ripristina lo Stack Pointer
		pop	ebp; ripristina il Base Pointer
		ret	   ; torna alla funzione C chiamante
