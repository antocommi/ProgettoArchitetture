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
extern norma1d

; ------------------------------------------------------------
; Funzione pagerank32
; ------------------------------------------------------------

global differenceNormd

;v1		equ		8       rdi 
;v2              equ             12     rsi  
;dim             equ             16     rdx
;dim2            equ             20     rcx
;delta           equ             24     r8

differenceNormd:
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq	                                        ; salva i registri generali

                xor r9, r9
                vxorpd ymm8, ymm8 
          ciclo:vmovapd ymm0, [rdi+r9*8]
                vmovapd ymm1, [rdi+r9*8+32]
                vmovapd ymm2, [rdi+r9*8+64]
                vmovapd ymm3, [rdi+r9*8+96]
                vmovapd ymm4, [rdi+r9*8+128]
                vmovapd ymm5, [rdi+r9*8+160]
                vmovapd ymm6, [rdi+r9*8+192]
                vmovapd ymm7, [rdi+r9*8+224]
                vmovapd ymm8, [rdi+r9*8+256]
                vmovapd ymm9, [rdi+r9*8+288]
                vmovapd ymm10, [rdi+r9*8+320]
                vmovapd ymm11, [rdi+r9*8+352]
                vmovapd ymm12, [rdi+r9*8+384]
                vmovapd ymm13, [rdi+r9*8+416]
                vmovapd ymm14, [rdi+r9*8+448]
                vmovapd ymm15, [rdi+r9*8+480]
                vsubpd ymm0, [rsi+r9*8]
                vsubpd ymm1, [rsi+r9*8+32] 
                vsubpd ymm2, [rsi+r9*8+64] 
                vsubpd ymm3, [rsi+r9*8+96] 
                vsubpd ymm4, [rsi+r9*8+128] 
                vsubpd ymm5, [rsi+r9*8+160] 
                vsubpd ymm6, [rsi+r9*8+192] 
                vsubpd ymm7, [rsi+r9*8+224]
                vsubpd ymm8, [rsi+r9*8+256]
                vsubpd ymm9, [rsi+r9*8+288] 
                vsubpd ymm10, [rsi+r9*8+320] 
                vsubpd ymm11, [rsi+r9*8+352] 
                vsubpd ymm12, [rsi+r9*8+384] 
                vsubpd ymm13, [rsi+r9*8+416] 
                vsubpd ymm14, [rsi+r9*8+448] 
                vsubpd ymm15, [rsi+r9*8+480]
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

                ;vhaddpd xmm8, xmm8   
                ;vperm2f128 ymm1, ymm8, ymm8, 1
                ;vaddsd xmm8, xmm1
                mov rsi, r8
                mov rdx, norma
                call norma1d
                vmovsd xmm8, [norma]
                vmovsd [rcx], xmm8
           

		popaq						; ripristina i registri generali
		mov		rsp, rbp			; ripristina lo Stack Pointer
		pop		rbp				; ripristina il Base Pointer
		ret						; torna alla funzione C chiamante
