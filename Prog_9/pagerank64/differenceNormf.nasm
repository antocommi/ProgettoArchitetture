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
; %include "norma1d.nasm"

section .data			; Sezione contenente dati inizializzati

uno:		dd		1.0
;
;align 16
;inizio:		dd		1.0, 2.0, 3.0, 4.0

section .bss			; Sezione contenente dati non inizializzati

alignb 16
;vec2:		resq	4
norma:          resd    1

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
extern norma1f

; ------------------------------------------------------------
; Funzione pagerank32
; ------------------------------------------------------------

global differenceNormf

;v1		equ		8       rdi 
;v2              equ             12     rsi  
;dim             equ             16     rdx
;dim2            equ             20     rcx
;delta           equ             24     r8

differenceNormf:
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq	                                        ; salva i registri generali

                xor r9, r9
                vxorps ymm8, ymm8
          ciclo:vmovaps ymm0, [rdi+r9*4]
                vmovaps ymm1, [rdi+r9*4+32]
                vmovaps ymm2, [rdi+r9*4+64]
                vmovaps ymm3, [rdi+r9*4+96]
                vmovaps ymm4, [rdi+r9*4+128]
                vmovaps ymm5, [rdi+r9*4+160]
                vmovaps ymm6, [rdi+r9*4+192]
                vmovaps ymm7, [rdi+r9*4+224]
                vmovaps ymm8, [rdi+r9*4+256]
                vmovaps ymm9, [rdi+r9*4+288]
                vmovaps ymm10, [rdi+r9*4+320]
                vmovaps ymm11, [rdi+r9*4+352]
                vmovaps ymm12, [rdi+r9*4+384]
                vmovaps ymm13, [rdi+r9*4+416]
                vmovaps ymm14, [rdi+r9*4+448]
                vmovaps ymm15, [rdi+r9*4+480]
                vsubps ymm0, [rsi+r9*4]
                vsubps ymm1, [rsi+r9*4+32] 
                vsubps ymm2, [rsi+r9*4+64] 
                vsubps ymm3, [rsi+r9*4+96] 
                vsubps ymm4, [rsi+r9*4+128] 
                vsubps ymm5, [rsi+r9*4+160] 
                vsubps ymm6, [rsi+r9*4+192] 
                vsubps ymm7, [rsi+r9*4+224]
                vsubps ymm8, [rsi+r9*4+256]
                vsubps ymm9, [rsi+r9*4+288] 
                vsubps ymm10, [rsi+r9*4+320] 
                vsubps ymm11, [rsi+r9*4+352] 
                vsubps ymm12, [rsi+r9*4+384] 
                vsubps ymm13, [rsi+r9*4+416] 
                vsubps ymm14, [rsi+r9*4+448] 
                vsubps ymm15, [rsi+r9*4+480]
                vmovaps [rdi+r9*4], ymm0
                vmovaps [rdi+r9*4+32], ymm1
                vmovaps [rdi+r9*4+64], ymm2
                vmovaps [rdi+r9*4+96], ymm3
                vmovaps [rdi+r9*4+128], ymm4
                vmovaps [rdi+r9*4+160], ymm5
                vmovaps [rdi+r9*4+192], ymm6
                vmovaps [rdi+r9*4+224], ymm7
                vmovaps [rdi+r9*4+256], ymm8
                vmovaps [rdi+r9*4+288], ymm9
                vmovaps [rdi+r9*4+320], ymm10
                vmovaps [rdi+r9*4+352], ymm11
                vmovaps [rdi+r9*4+384], ymm12
                vmovaps [rdi+r9*4+416], ymm13
                vmovaps [rdi+r9*4+448], ymm14
                vmovaps [rdi+r9*4+480], ymm15
                add r9, 128
                cmp r9, rdx
                jl ciclo

                ;vhaddps xmm8, xmm8  
                ;vhaddps xmm8, xmm8 
                ;vperm2f128 ymm1, ymm8, ymm8, 1
                ;vaddss xmm8, xmm1
                mov rsi, r8
                mov rdx, norma
                call norma1f
                vmovss xmm8, [norma]
                vmovss [rcx], xmm8
           

		popaq						; ripristina i registri generali
		mov		rsp, rbp			; ripristina lo Stack Pointer
		pop		rbp				; ripristina il Base Pointer
		ret						; torna alla funzione C chiamante
