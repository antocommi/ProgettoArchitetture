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

;alignb 16
;vec2:		resq	4

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

global assegnaD

;v1		equ		8     rdi
;v2              equ             12   rsi
;n               equ             16   rdx

assegnaD:
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq	                                        ; salva i registri generali

    
                xor r9, r9
          ciclo:vmovapd ymm0, [rsi+r9*8]
                vmovapd ymm1, [rsi+r9*8+32]
                vmovapd ymm2, [rsi+r9*8+64]
                vmovapd ymm3, [rsi+r9*8+96]
                vmovapd ymm4, [rsi+r9*8+128]
                vmovapd ymm5, [rsi+r9*8+160]
                vmovapd ymm6, [rsi+r9*8+192]
                vmovapd ymm7, [rsi+r9*8+224]
                vmovapd ymm8, [rsi+r9*8+256]
                vmovapd ymm9, [rsi+r9*8+288]
                vmovapd ymm10, [rsi+r9*8+320]
                vmovapd ymm11, [rsi+r9*8+352]
                vmovapd ymm12, [rsi+r9*8+384]
                vmovapd ymm13, [rsi+r9*8+416]
                vmovapd ymm14, [rsi+r9*8+448]
                vmovapd ymm15, [rsi+r9*8+480]
                vmovapd [rdi+r9*8], ymm0
                vmovapd [rdi+r9*8+32], ymm1
                vmovapd [rdi+r9*8+64], ymm2
                vmovapd [rdi+r9*8+96], ymm3
                vmovapd [rdi+r9*8+128], ymm4
                vmovapd [rdi+r9*8+160], ymm5
                vmovapd [rdi+r9*8+192], ymm6
                vmovapd [rdi+r9*8+224], ymm7
                vmovapd [rdi+r9*8+256], ymm8
                vmovapd [rdi+r9*8+288], ymm9
                vmovapd [rdi+r9*8+320], ymm10
                vmovapd [rdi+r9*8+352], ymm11
                vmovapd [rdi+r9*8+384], ymm12
                vmovapd [rdi+r9*8+416], ymm13
                vmovapd [rdi+r9*8+448], ymm14
                vmovapd [rdi+r9*8+480], ymm15
                add r9, 64
                cmp r9, rdx
                jl ciclo

		popaq						; ripristina i registri generali
		mov		rsp, rbp			; ripristina lo Stack Pointer
		pop		rbp				; ripristina il Base Pointer
		ret						; torna alla funzione C chiamante
