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

;l		equ		8   rdi
;temp            equ             12 rsi
;p               equ             16 xmm0
;n               equ             24 rdx

powerExtrd:
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq
		
		
                xor r9, r9
                vinsertf128 ymm0, ymm0, xmm0, 1 
                ;vpermilpd ymm0, ymm0, 0
                vshufpd ymm0, ymm0, ymm0, 0
                vbroadcastsd ymm1, [uno]
                vsubpd ymm1, ymm0
          ciclo:vmovapd ymm2, [rdi+r9*8]
                vmovapd ymm3, [rdi+r9*8+32]
                vmovapd ymm4, [rdi+r9*8+64]
                vmovapd ymm5, [rdi+r9*8+96]
                vmulpd ymm2, ymm0
                vmulpd ymm3, ymm0
                vmulpd ymm4, ymm0
                vmulpd ymm5, ymm0
                vmovapd ymm6, [rsi+r9*8]
                vmovapd ymm7, [rsi+r9*8+32]
                vmovapd ymm8, [rsi+r9*8+64]
                vmovapd ymm9, [rsi+r9*8+96]
                vsubpd ymm6, ymm2
                vsubpd ymm7, ymm3
                vsubpd ymm8, ymm4
                vsubpd ymm9, ymm5
                vdivpd ymm6, ymm1
                vdivpd ymm7, ymm1
                vdivpd ymm8, ymm1
                vdivpd ymm9, ymm1
                vmovapd [rsi+r9*8], ymm6
                vmovapd [rsi+r9*8+32], ymm7
                vmovapd [rsi+r9*8+64], ymm8
                vmovapd [rsi+r9*8+96], ymm9
                vmovapd ymm2, [rdi+r9*8+128]
                vmovapd ymm3, [rdi+r9*8+160]
                vmovapd ymm4, [rdi+r9*8+192]
                vmovapd ymm5, [rdi+r9*8+224]
                vmulpd ymm2, ymm0
                vmulpd ymm3, ymm0
                vmulpd ymm4, ymm0
                vmulpd ymm5, ymm0
                vmovapd ymm6, [rsi+r9*8+128]
                vmovapd ymm7, [rsi+r9*8+160]
                vmovapd ymm8, [rsi+r9*8+192]
                vmovapd ymm9, [rsi+r9*8+224]
                vsubpd ymm6, ymm2
                vsubpd ymm7, ymm3
                vsubpd ymm8, ymm4
                vsubpd ymm9, ymm5
                vdivpd ymm6, ymm1
                vdivpd ymm7, ymm1
                vdivpd ymm8, ymm1
                vdivpd ymm9, ymm1
                vmovapd [rsi+r9*8+128], ymm6
                vmovapd [rsi+r9*8+160], ymm7
                vmovapd [rsi+r9*8+192], ymm8
                vmovapd [rsi+r9*8+224], ymm9
                add r9, 32
                cmp r9, rdx
                jl ciclo

		popaq						; ripristina i registri generali
		mov		rsp, rbp			; ripristina lo Stack Pointer
		pop		rbp				; ripristina il Base Pointer
		ret						; torna alla funzione C chiamante
