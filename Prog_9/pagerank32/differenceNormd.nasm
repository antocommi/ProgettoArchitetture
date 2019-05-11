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
; %include "norma1d.nasm"

section .data			; Sezione contenente dati inizializzati

uno:		dd		1.0
;
;align 16
;inizio:		dd		1.0, 2.0, 3.0, 4.0

section .bss			; Sezione contenente dati non inizializzati

alignb 16
;vec2:		resq	4
norma:          resq    1

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

; ------------------------------------------------------------
; Funzione pagerank32
; ------------------------------------------------------------

global differenceNormd

v1		equ		8
v2              equ             12
dim             equ             16
delta           equ             20

differenceNormd:
		push		ebp; salva il Base Pointer
		mov		ebp, esp; il Base Pointer punta al Record di Attivazione corrente
		push		ebx; salva i registri da preservare
		push		esi
		push		edi
		
		mov		eax, [ebp+v1]
                mov             ebx, [ebp+v2]
                mov             edx, [ebp+dim]
                mov             ecx, [ebp+delta]
  
                xor esi, esi
                xorpd xmm4, xmm4
          ciclo:movapd xmm0, [eax+esi*8]
                movapd xmm1, [eax+esi*8+16]
                movapd xmm2, [eax+esi*8+32]
                movapd xmm3, [eax+esi*8+48]
                subpd xmm0, [ebx+esi*8]
                subpd xmm1, [ebx+esi*8+16] 
                subpd xmm2, [ebx+esi*8+32] 
                subpd xmm3, [ebx+esi*8+48] 
                psllq xmm0,1
                psllq xmm1,1
                psllq xmm2,1
                psllq xmm3,1
                psrlq xmm0,1
                psrlq xmm1,1
                psrlq xmm2,1
                psrlq xmm3,1
                addpd xmm4,xmm1
                addpd xmm4,xmm2
                addpd xmm4,xmm3
                addpd xmm4,xmm0
                add esi,8
                cmp esi,edx
                jl ciclo

               
                haddpd xmm4, xmm4
                ;movhlps xmm1, xmm4
                ;addsd xmm4, xmm1
                movsd [ecx], xmm4

		pop	edi; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp; ripristina lo Stack Pointer
		pop	ebp; ripristina il Base Pointer
		ret	   ; torna alla funzione C chiamante
