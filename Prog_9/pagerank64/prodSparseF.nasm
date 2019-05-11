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
;align 16
uno:		dq		1.0
;
;align 16
;inizio:		dd		1.0, 2.0, 3.0, 4.0

section .bss			; Sezione contenente dati non inizializzati

alignb 32
somma:		resq	1
rii:             resq    4

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
extern sommad

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

global prodSparseF

;row		equ		8    rdi
;col             equ             12  rsi
;vd              equ             16  rdx
;N               equ             20  rcx
;xk1             equ             24  r8
;xk              equ             28  r9
d               equ             16
M               equ             24
prodSparseF:
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq	                                        ; salva i registri generali
  
                mov rax, [rbp+d]
                mov rbx, [rbp+M] 

                xor r10, r10
         cicloV:vmovaps ymm0,[r8+r10*4]
                vmovaps ymm1,[r8+r10*4+32]
                vmovaps ymm2,[r8+r10*4+64]
                vmovaps ymm3,[r8+r10*4+96]
                vmovaps ymm4,[r8+r10*4+128]
                vmovaps ymm5,[r8+r10*4+160]
                vmovaps ymm6,[r8+r10*4+192]
                vmovaps ymm7,[r8+r10*4+224]
                vmovaps ymm8,[r8+r10*4+256]
                vmovaps ymm9,[r8+r10*4+288]
                vmovaps ymm10,[r8+r10*4+320]
                vmovaps ymm11,[r8+r10*4+352]
                vmovaps ymm12,[r8+r10*4+384]
                vmovaps ymm13,[r8+r10*4+416]
                vmovaps ymm14,[r8+r10*4+448]
                vmovaps ymm15,[r8+r10*4+480]
                vmulps ymm0,[rax+r10*4]
                vmulps ymm1,[rax+r10*4+32]
                vmulps ymm2,[rax+r10*4+64]
                vmulps ymm3,[rax+r10*4+96]
                vmulps ymm4,[rax+r10*4+128]
                vmulps ymm5,[rax+r10*4+160]
                vmulps ymm6,[rax+r10*4+192]
                vmulps ymm7,[rax+r10*4+224]
                vmulps ymm8,[rax+r10*4+256]
                vmulps ymm9,[rax+r10*4+288]
                vmulps ymm10,[rax+r10*4+320]
                vmulps ymm11,[rax+r10*4+352]
                vmulps ymm12,[rax+r10*4+384]
                vmulps ymm13,[rax+r10*4+416]
                vmulps ymm14,[rax+r10*4+448]
                vmulps ymm15,[rax+r10*4+480]
                vmovaps [rdx+r10*4],ymm0
                vmovaps [rdx+r10*4+32],ymm1
                vmovaps [rdx+r10*4+64],ymm2
                vmovaps [rdx+r10*4+96],ymm3
                vmovaps [rdx+r10*4+128],ymm4
                vmovaps [rdx+r10*4+160],ymm5
                vmovaps [rdx+r10*4+192],ymm6
                vmovaps [rdx+r10*4+224],ymm7
                vmovaps [rdx+r10*4+256],ymm8
                vmovaps [rdx+r10*4+288],ymm9
                vmovaps [rdx+r10*4+320],ymm10
                vmovaps [rdx+r10*4+352],ymm11
                vmovaps [rdx+r10*4+384],ymm12
                vmovaps [rdx+r10*4+416],ymm13
                vmovaps [rdx+r10*4+448],ymm14
                vmovaps [rdx+r10*4+480],ymm15
                add r10, 128
                cmp r10, rbx
                jl cicloV
                
                xor r10, r10
                xor r11, r11
                mov r12d, [rsi]
          ciclo:
                vxorps ymm0, ymm0
                cmp r12, 0
                je fine
                
         ciclor:
                mov r13d, [rdi+r11*4]
                sal r13d, 2
                add r13, rdx
                vmovss xmm1, [r13]
                vaddss xmm0, xmm1
   
                inc r11
                dec r12
                cmp r12, 0
                jg ciclor
                
           fine:
                vmovss [r9+r10*4], xmm0
                mov r12d, [rsi+r10*4+4]
                sub r12d, [rsi+r10*4]
                inc r10  
                cmp r10, rcx
                jl ciclo
                

		popaq						; ripristina i registri generali
		mov		rsp, rbp			; ripristina lo Stack Pointer
		pop		rbp				; ripristina il Base Pointer
		ret						; torna alla funzione C chiamante
