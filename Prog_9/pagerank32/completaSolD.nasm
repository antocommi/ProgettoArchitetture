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

uno:		dd		1.0
;
;align 16
;inizio:		dd		1.0, 2.0, 3.0, 4.0

section .bss			; Sezione contenente dati non inizializzati

;alignb 16
;vec2:		resq	4
ri: resq 1

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

global completaSolD

xk		equ		8
dim             equ             12
val             equ             16

completaSolD:
		push		ebp; salva il Base Pointer
		mov		ebp, esp; il Base Pointer punta al Record di Attivazione corrente
		push		ebx; salva i registri da preservare
		push		esi
		push		edi
		
		mov		eax, [ebp+dim]
                mov             ebx, [ebp+xk]
              
                xor esi, esi
                movsd xmm0, [ebp+val]
                shufpd xmm0, xmm0, 0
         ciclo2:movapd xmm1, xmm0
                movapd xmm2, xmm0
                movapd xmm3, xmm0
                movapd xmm4, xmm0
                addpd xmm1, [ebx+esi*8]
                addpd xmm2, [ebx+esi*8+16]
                addpd xmm3, [ebx+esi*8+32]
                addpd xmm4, [ebx+esi*8+48]
                movapd [ebx+esi*8], xmm1
                movapd [ebx+esi*8+16], xmm2
                movapd [ebx+esi*8+32], xmm3
                movapd [ebx+esi*8+48], xmm4
                movapd xmm1, xmm0
                movapd xmm2, xmm0
                movapd xmm3, xmm0
                movapd xmm4, xmm0
                addpd xmm1, [ebx+esi*8+64]
                addpd xmm2, [ebx+esi*8+80]
                addpd xmm3, [ebx+esi*8+96]
                addpd xmm4, [ebx+esi*8+112]
                movapd [ebx+esi*8+64], xmm1
                movapd [ebx+esi*8+80], xmm2
                movapd [ebx+esi*8+96], xmm3
                movapd [ebx+esi*8+112], xmm4
                add esi, 16
                cmp esi, eax
                jl ciclo2

		pop	edi; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp; ripristina lo Stack Pointer
		pop	ebp; ripristina il Base Pointer
		ret	   ; torna alla funzione C chiamante
