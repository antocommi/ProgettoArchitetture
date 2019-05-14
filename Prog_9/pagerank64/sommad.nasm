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

global sommad

;v		equ		8     rdi
;dim             equ             12   rsi
;s               equ             16   rdx

sommad:
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq	                                        ; salva i registri generali


                vmovapd ymm0, [rdi]
                vmovapd ymm1, [rdi+32]
                vmovapd ymm2, [rdi+64]
                vmovapd ymm3, [rdi+96]
                vmovapd ymm4, [rdi+128]
                vmovapd ymm5, [rdi+160]
                vmovapd ymm6, [rdi+192]
                vmovapd ymm7, [rdi+224]
                vmovapd ymm8, [rdi+256]
                vmovapd ymm9, [rdi+288]
                vmovapd ymm10, [rdi+320]
                vmovapd ymm11, [rdi+352]
                vmovapd ymm12, [rdi+384]
                vmovapd ymm13, [rdi+416]
                vmovapd ymm14, [rdi+448]
                vmovapd ymm15, [rdi+480]
                cmp rsi, 64
                je fine
                mov r8, 64
         cicloq:vaddpd ymm0, [rdi+r8*8]
                vaddpd ymm1, [rdi+r8*8+32]
                vaddpd ymm2, [rdi+r8*8+64]
                vaddpd ymm3, [rdi+r8*8+96]
                vaddpd ymm4, [rdi+r8*8+128]
                vaddpd ymm5, [rdi+r8*8+160]
                vaddpd ymm6, [rdi+r8*8+192]
                vaddpd ymm7, [rdi+r8*8+224]
                vaddpd ymm8, [rdi+r8*8+256]
                vaddpd ymm9, [rdi+r8*8+288]
                vaddpd ymm10, [rdi+r8*8+320]
                vaddpd ymm11, [rdi+r8*8+352]
                vaddpd ymm12, [rdi+r8*8+384]
                vaddpd ymm13, [rdi+r8*8+416]
                vaddpd ymm14, [rdi+r8*8+448]
                vaddpd ymm15, [rdi+r8*8+480]
                add r8, 64
                cmp r8, rsi
                jl cicloq 
           fine:vaddpd ymm0, ymm1
                vaddpd ymm0, ymm2
                vaddpd ymm0, ymm3
                vaddpd ymm0, ymm4
                vaddpd ymm0, ymm5
                vaddpd ymm0, ymm6
                vaddpd ymm0, ymm7
                vaddpd ymm0, ymm8
                vaddpd ymm0, ymm9
                vaddpd ymm0, ymm10
                vaddpd ymm0, ymm11
                vaddpd ymm0, ymm12
                vaddpd ymm0, ymm13
                vaddpd ymm0, ymm14
                vaddpd ymm0, ymm15
                vhaddpd ymm0, ymm0
                vperm2f128 ymm1, ymm0, ymm0, 1
                vaddsd xmm0, xmm1
                movsd [rdx], xmm0	
                

		popaq						; ripristina i registri generali
		mov		rsp, rbp			; ripristina lo Stack Pointer
		pop		rbp				; ripristina il Base Pointer
		ret						; torna alla funzione C chiamante
