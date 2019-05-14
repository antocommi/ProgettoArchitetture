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
; %include "norma1f.nasm"

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

global differenceNormf

v1		equ		8
v2              equ             12
dim             equ             16
delta           equ             20

differenceNormf:
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
                xorps xmm4, xmm4
          ciclo:movaps xmm0, [eax+esi*4]
                movaps xmm1, [eax+esi*4+16]
                movaps xmm2, [eax+esi*4+32]
                movaps xmm3, [eax+esi*4+48]
                subps xmm0, [ebx+esi*4]
                subps xmm1, [ebx+esi*4+16] 
                subps xmm2, [ebx+esi*4+32] 
                subps xmm3, [ebx+esi*4+48] 
                pslld xmm0,1
                pslld xmm1,1
                pslld xmm2,1
                pslld xmm3,1
                psrld xmm0,1
                psrld xmm1,1
                psrld xmm2,1
                psrld xmm3,1
                addps xmm4,xmm1
                addps xmm4,xmm2
                addps xmm4,xmm3
                addps xmm4,xmm0
                add esi,16
                cmp esi,edx
                jl ciclo

                haddps xmm4, xmm4   
                haddps xmm4, xmm4
                movss [ecx], xmm4
           

		pop	edi; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp; ripristina lo Stack Pointer
		pop	ebp; ripristina il Base Pointer
		ret	   ; torna alla funzione C chiamante
