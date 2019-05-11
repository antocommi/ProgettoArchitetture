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
R:		resd	1

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

global sommad

v		equ		8
dim             equ             12
s               equ             16

sommad:
		push		ebp; salva il Base Pointer
		mov		ebp, esp; il Base Pointer punta al Record di Attivazione corrente
		push		ebx; salva i registri da preservare
		push		esi
		push		edi
		
		mov		eax, [ebp+v]
                mov             ebx, [ebp+dim]
                mov             ecx, [ebp+s]

                movapd xmm0, [eax]
                movapd xmm1, [eax+16]
                movapd xmm2, [eax+32]
                movapd xmm3, [eax+48]
                movapd xmm4, [eax+64]
                movapd xmm5, [eax+80]
                movapd xmm6, [eax+96]
                movapd xmm7, [eax+112]
                cmp ebx,16
                je fine
                mov esi,16
         cicloq:addpd xmm0, [eax+esi*8]
                addpd xmm1, [eax+esi*8+16]
                addpd xmm2, [eax+esi*8+32]
                addpd xmm3, [eax+esi*8+48]
                addpd xmm4, [eax+esi*8+64]
                addpd xmm5, [eax+esi*8+80]
                addpd xmm6, [eax+esi*8+96]
                addpd xmm7, [eax+esi*8+112]
                add esi,16
                cmp esi,ebx
                jl cicloq 
           fine:addpd xmm0, xmm1
                addpd xmm0, xmm2
                addpd xmm0, xmm3
                addpd xmm0, xmm4
                addpd xmm0, xmm5
                addpd xmm0, xmm6
                addpd xmm0, xmm7
                haddpd xmm0, xmm0
                movsd [ecx], xmm0	
                

		pop	edi; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp; ripristina lo Stack Pointer
		pop	ebp; ripristina il Base Pointer
		ret	   ; torna alla funzione C chiamante
