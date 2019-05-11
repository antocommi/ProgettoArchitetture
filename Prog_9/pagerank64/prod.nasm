; ---------------------------------------------------------
; PageRank con istruzioni AVX a 64 bit
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
;     nasm -f elf64 pagerank64.nasm 
;

%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati

zero:		dq		0.0
;
;align 32
;vec1:		dd		1.0, 2.0, 3.0, 4.0

section .bss			; Sezione contenente dati non inizializzati

;alignb 32
;vec2:		resq	4
ri:    resq 1
ri2:  resq 4

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
	mov	rdi, %1
	mov	rsi, %2
	call	get_block
%endmacro

%macro	fremem	1
	mov	rdi, %1
	call	free_block
%endmacro



; ------------------------------------------------------------
; Funzione pagerank64
; ------------------------------------------------------------

global prod

;P		equ		8     rdi
;v               equ             12   rsi
;dim             equ             16   rdx
;n               equ             20   rcx
;temp            equ             24   r8

prod:

		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq						; salva i registri generali


		xor r9, r9
                xor r10, r10
         cicloi:xor r11, r11
                vxorpd ymm0, ymm0
                vxorpd ymm1, ymm1
                vxorpd ymm2, ymm2
                vxorpd ymm3, ymm3
                vxorpd ymm8, ymm8
                vxorpd ymm9, ymm9
                vxorpd ymm10, ymm10
                vxorpd ymm11, ymm11
         cicloj:vmovapd ymm4, [rdi+r10*8]
                vmovapd ymm5, [rdi+r10*8+32] 
                vmovapd ymm6, [rdi+r10*8+64] 
                vmovapd ymm7, [rdi+r10*8+96]
                vmovapd ymm12, [rdi+r10*8+128]
                vmovapd ymm13, [rdi+r10*8+160] 
                vmovapd ymm14, [rdi+r10*8+192] 
                vmovapd ymm15, [rdi+r10*8+224]
                vmulpd ymm4, [rsi+r11*8]
                vmulpd ymm5, [rsi+r11*8+32]
                vmulpd ymm6, [rsi+r11*8+64]
                vmulpd ymm7, [rsi+r11*8+96]
                vmulpd ymm12, [rsi+r11*8+128]
                vmulpd ymm13, [rsi+r11*8+160] 
                vmulpd ymm14, [rsi+r11*8+192] 
                vmulpd ymm15, [rsi+r11*8+224]
                vaddpd ymm0, ymm4
                vaddpd ymm1, ymm5
                vaddpd ymm2, ymm6
                vaddpd ymm3, ymm7
                vaddpd ymm8, ymm12
                vaddpd ymm9, ymm13
                vaddpd ymm10, ymm14
                vaddpd ymm11, ymm15
                add r10, 32
                add r11, 32
                cmp r11, rdx
                jl cicloj
                vaddpd ymm0, ymm1
                vaddpd ymm0, ymm2
                vaddpd ymm0, ymm3
                vaddpd ymm0, ymm8
                vaddpd ymm0, ymm9
                vaddpd ymm0, ymm10
                vaddpd ymm0, ymm11
                vhaddpd ymm0, ymm0
                vperm2f128 ymm1, ymm0, ymm0, 1
                vaddsd xmm0, xmm1
                vmovsd [r8+r9*8], xmm0
                inc r9
                cmp r9, rcx
                jl cicloi  
		
		
		
		popaq						; ripristina i registri generali
		mov		rsp, rbp			; ripristina lo Stack Pointer
		pop		rbp				; ripristina il Base Pointer
		ret						; torna alla funzione C chiamante
		
