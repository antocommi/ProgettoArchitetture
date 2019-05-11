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

uno:		dq		1.0
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

global powerExtrd

l		equ		8
temp            equ             12
p               equ             16
n               equ             24

powerExtrd:
		push		ebp; salva il Base Pointer
		mov		ebp, esp; il Base Pointer punta al Record di Attivazione corrente
		push		ebx; salva i registri da preservare
		push		esi
		push		edi
		
		mov		eax, [ebp+l]
                mov             ebx, [ebp+temp]
                mov             edi, [ebp+n]

                xor esi, esi
                movsd xmm6, [ebp+p]
                shufpd xmm6, xmm6, 0
                movsd xmm7, [uno]
                subsd xmm7, xmm6
                shufpd xmm7, xmm7, 0
          ciclo:movapd xmm0, [eax+esi*8]
                movapd xmm1, [eax+esi*8+16]
                mulpd xmm0, xmm6
                mulpd xmm1, xmm6
                movapd xmm2, [ebx+esi*8]
                movapd xmm3, [ebx+esi*8+16]
                subpd xmm2, xmm0
                subpd xmm3, xmm1
                divpd xmm2, xmm7
                divpd xmm3, xmm7
                movapd [ebx+esi*8], xmm2
                movapd [ebx+esi*8+16], xmm3
                movapd xmm0, [eax+esi*8+32]
                movapd xmm1, [eax+esi*8+48]
                mulpd xmm0, xmm6
                mulpd xmm1, xmm6
                movapd xmm2, [ebx+esi*8+32]
                movapd xmm3, [ebx+esi*8+48]
                subpd xmm2, xmm0
                subpd xmm3, xmm1
                divpd xmm2, xmm7
                divpd xmm3, xmm7
                movapd [ebx+esi*8+32], xmm2
                movapd [ebx+esi*8+48], xmm3              
                add esi, 8
                cmp esi, edi
                jl ciclo

		pop	edi; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp; ripristina lo Stack Pointer
		pop	ebp; ripristina il Base Pointer
		ret	   ; torna alla funzione C chiamante
