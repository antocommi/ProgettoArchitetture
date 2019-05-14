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

 %include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati

uno:		dd		1.0
;
;align 16
;inizio:		dd		1.0, 2.0, 3.0, 4.0

section .bss			; Sezione contenente dati non inizializzati

alignb 32
;vec2:		resq	4
ri: resq 4

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

global completaSolF

;xk		equ		8   rdi
;dim             equ             12 rsi
;val             equ             16 xmm0

completaSolF:
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq	                                        ; salva i registri generali

              
                xor r9, r9
                vinsertf128 ymm0, ymm0, xmm0, 1 
                vpermilps ymm0, ymm0, 0
         ciclo2:vmovaps ymm1, ymm0
                vmovaps ymm2, ymm0
                vmovaps ymm3, ymm0
                vmovaps ymm4, ymm0
                vmovaps ymm5, ymm0
                vmovaps ymm6, ymm0
                vmovaps ymm7, ymm0
                vmovaps ymm8, ymm0
                vaddps ymm1, [rdi+r9*4]
                vaddps ymm2, [rdi+r9*4+32]
                vaddps ymm3, [rdi+r9*4+64]
                vaddps ymm4, [rdi+r9*4+96]
                vaddps ymm5, [rdi+r9*4+128]
                vaddps ymm6, [rdi+r9*4+160]
                vaddps ymm7, [rdi+r9*4+192]
                vaddps ymm8, [rdi+r9*4+224]
                vmovaps [rdi+r9*4], ymm1
                vmovaps [rdi+r9*4+32], ymm2
                vmovaps [rdi+r9*4+64], ymm3
                vmovaps [rdi+r9*4+96], ymm4
                vmovaps [rdi+r9*4+128], ymm5
                vmovaps [rdi+r9*4+160], ymm6
                vmovaps [rdi+r9*4+192], ymm7
                vmovaps [rdi+r9*4+224], ymm8
                vmovaps ymm1, ymm0
                vmovaps ymm2, ymm0
                vmovaps ymm3, ymm0
                vmovaps ymm4, ymm0
                vmovaps ymm5, ymm0
                vmovaps ymm6, ymm0
                vmovaps ymm7, ymm0
                vmovaps ymm8, ymm0
                vaddps ymm1, [rdi+r9*4+256]
                vaddps ymm2, [rdi+r9*4+288]
                vaddps ymm3, [rdi+r9*4+320]
                vaddps ymm4, [rdi+r9*4+352]
                vaddps ymm5, [rdi+r9*4+384]
                vaddps ymm6, [rdi+r9*4+416]
                vaddps ymm7, [rdi+r9*4+448]
                vaddps ymm8, [rdi+r9*4+480]
                vmovaps [rdi+r9*4+256], ymm1
                vmovaps [rdi+r9*4+288], ymm2
                vmovaps [rdi+r9*4+320], ymm3
                vmovaps [rdi+r9*4+352], ymm4
                vmovaps [rdi+r9*4+384], ymm5
                vmovaps [rdi+r9*4+416], ymm6
                vmovaps [rdi+r9*4+448], ymm7
                vmovaps [rdi+r9*4+480], ymm8
                add r9, 128
                cmp r9, rsi
                jl ciclo2

		popaq						; ripristina i registri generali
		mov		rsp, rbp			; ripristina lo Stack Pointer
		pop		rbp				; ripristina il Base Pointer
		ret						; torna alla funzione C chiamante
