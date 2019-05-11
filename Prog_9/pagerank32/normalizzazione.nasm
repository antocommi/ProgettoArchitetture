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
;%include "norma1d.nasm"

section .data			; Sezione contenente dati inizializzati

uno:		dd		1.0
;
;align 16
;inizio:		dd		1.0, 2.0, 3.0, 4.0

section .bss			; Sezione contenente dati non inizializzati

;alignb 16
norma:		resq	1
ri: resd 1

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

global normalizzazione

v		equ		8
dim             equ             12
dim2            equ             16

normalizzazione:
		push		ebp; salva il Base Pointer
		mov		ebp, esp; il Base Pointer punta al Record di Attivazione corrente
		push		ebx; salva i registri da preservare
		push		esi
		push		edi
		
		mov		eax, [ebp+v]
                mov             ebx, [ebp+dim]
                mov             edi, [ebp+dim2]
                
                ;pushad
                ;push norma
                ;push ebx
                ;push eax
                ;call norma1d
                ;add esp, 12
                ;popad
               xor esi, esi
               xorpd xmm0, xmm0
        ciclo1:
               movapd xmm1, [eax+esi*8]
               movapd xmm2, [eax+esi*8+16]
               movapd xmm3, [eax+esi*8+32]
               movapd xmm4, [eax+esi*8+48]
               psllq xmm4,1
               psllq xmm1,1
               psllq xmm2,1
               psllq xmm3,1
               psrlq xmm4,1
               psrlq xmm1,1
               psrlq xmm2,1
               psrlq xmm3,1
               addpd xmm0,xmm1
               addpd xmm0,xmm2
               addpd xmm0,xmm3
               addpd xmm0,xmm4
               add esi,8
               cmp esi,ebx
               jl ciclo1

               
               haddpd xmm0, xmm0
                shufpd xmm0, xmm0, 0
                xor esi, esi
          ciclo:movapd xmm1, [eax+esi*8]
                movapd xmm2, [eax+esi*8+16]
                movapd xmm3, [eax+esi*8+32]
                movapd xmm4, [eax+esi*8+48]
                divpd xmm1, xmm0
                divpd xmm2, xmm0
                divpd xmm3, xmm0
                divpd xmm4, xmm0
                movapd [eax+esi*8], xmm1
                movapd [eax+esi*8+16], xmm2
                movapd [eax+esi*8+32], xmm3
                movapd [eax+esi*8+48], xmm4
                movapd xmm1, [eax+esi*8+64]
                movapd xmm2, [eax+esi*8+80]
                movapd xmm3, [eax+esi*8+96]
                movapd xmm4, [eax+esi*8+112]
                divpd xmm1, xmm0
                divpd xmm2, xmm0
                divpd xmm3, xmm0
                divpd xmm4, xmm0
                movapd [eax+esi*8+64], xmm1
                movapd [eax+esi*8+80], xmm2
                movapd [eax+esi*8+96], xmm3
                movapd [eax+esi*8+112], xmm4
                add esi, 16
                cmp esi, edi
                jl ciclo

		pop	edi; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp; ripristina lo Stack Pointer
		pop	ebp; ripristina il Base Pointer
		ret	   ; torna alla funzione C chiamante
