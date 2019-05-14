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

global sommaf

;v		equ		8     rdi
;dim             equ             12   rsi
;s               equ             16   rdx

sommaf:
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq	                                        ; salva i registri generali


                vmovaps ymm0, [rdi]
                vmovaps ymm1, [rdi+32]
                vmovaps ymm2, [rdi+64]
                vmovaps ymm3, [rdi+96]
                vmovaps ymm4, [rdi+128]
                vmovaps ymm5, [rdi+160]
                vmovaps ymm6, [rdi+192]
                vmovaps ymm7, [rdi+224]
                vmovaps ymm8, [rdi+256]
                vmovaps ymm9, [rdi+288]
                vmovaps ymm10, [rdi+320]
                vmovaps ymm11, [rdi+352]
                vmovaps ymm12, [rdi+384]
                vmovaps ymm13, [rdi+416]
                vmovaps ymm14, [rdi+448]
                vmovaps ymm15, [rdi+480]
                cmp rsi, 128
                je fine
                mov r8, 128
         cicloq:vaddps ymm0, [rdi+r8*4]
                vaddps ymm1, [rdi+r8*4+32]
                vaddps ymm2, [rdi+r8*4+64]
                vaddps ymm3, [rdi+r8*4+96]
                vaddps ymm4, [rdi+r8*4+128]
                vaddps ymm5, [rdi+r8*4+160]
                vaddps ymm6, [rdi+r8*4+192]
                vaddps ymm7, [rdi+r8*4+224]
                vaddps ymm8, [rdi+r8*4+256]
                vaddps ymm9, [rdi+r8*4+288]
                vaddps ymm10, [rdi+r8*4+320]
                vaddps ymm11, [rdi+r8*4+352]
                vaddps ymm12, [rdi+r8*4+384]
                vaddps ymm13, [rdi+r8*4+416]
                vaddps ymm14, [rdi+r8*4+448]
                vaddps ymm15, [rdi+r8*4+480]
                add r8, 128
                cmp r8, rsi
                jl cicloq 
           fine:vaddps ymm0, ymm1
                vaddps ymm0, ymm2
                vaddps ymm0, ymm3
                vaddps ymm0, ymm4
                vaddps ymm0, ymm5
                vaddps ymm0, ymm6
                vaddps ymm0, ymm7
                vaddps ymm0, ymm8
                vaddps ymm0, ymm9
                vaddps ymm0, ymm10
                vaddps ymm0, ymm11
                vaddps ymm0, ymm12
                vaddps ymm0, ymm13
                vaddps ymm0, ymm14
                vaddps ymm0, ymm15
                vhaddps ymm0, ymm0
                vhaddps ymm0, ymm0
                vperm2f128 ymm1, ymm0, ymm0, 1
                vaddss xmm0, xmm1
                movss [rdx], xmm0	
                

		popaq						; ripristina i registri generali
		mov		rsp, rbp			; ripristina lo Stack Pointer
		pop		rbp				; ripristina il Base Pointer
		ret						; torna alla funzione C chiamante
