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

global prodSparseD

;row		equ		8    rdi
;col             equ             12  rsi
;vd              equ             16  rdx
;N               equ             20  rcx
;xk1             equ             24  r8
;xk              equ             28  r9
d               equ             16
M               equ             24
prodSparseD:
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq	                                        ; salva i registri generali
  
                mov rax, [rbp+d]
                mov rbx, [rbp+M] 

                xor r10, r10
         cicloV:vmovapd ymm0,[r8+r10*8]
                vmovapd ymm1,[r8+r10*8+32]
                vmovapd ymm2,[r8+r10*8+64]
                vmovapd ymm3,[r8+r10*8+96]
                vmovapd ymm4,[r8+r10*8+128]
                vmovapd ymm5,[r8+r10*8+160]
                vmovapd ymm6,[r8+r10*8+192]
                vmovapd ymm7,[r8+r10*8+224]
                vmovapd ymm8,[r8+r10*8+256]
                vmovapd ymm9,[r8+r10*8+288]
                vmovapd ymm10,[r8+r10*8+320]
                vmovapd ymm11,[r8+r10*8+352]
                vmovapd ymm12,[r8+r10*8+384]
                vmovapd ymm13,[r8+r10*8+416]
                vmovapd ymm14,[r8+r10*8+448]
                vmovapd ymm15,[r8+r10*8+480]
                vmulpd ymm0,[rax+r10*8]
                vmulpd ymm1,[rax+r10*8+32]
                vmulpd ymm2,[rax+r10*8+64]
                vmulpd ymm3,[rax+r10*8+96]
                vmulpd ymm4,[rax+r10*8+128]
                vmulpd ymm5,[rax+r10*8+160]
                vmulpd ymm6,[rax+r10*8+192]
                vmulpd ymm7,[rax+r10*8+224]
                vmulpd ymm8,[rax+r10*8+256]
                vmulpd ymm9,[rax+r10*8+288]
                vmulpd ymm10,[rax+r10*8+320]
                vmulpd ymm11,[rax+r10*8+352]
                vmulpd ymm12,[rax+r10*8+384]
                vmulpd ymm13,[rax+r10*8+416]
                vmulpd ymm14,[rax+r10*8+448]
                vmulpd ymm15,[rax+r10*8+480]
                vmovapd [rdx+r10*8],ymm0
                vmovapd [rdx+r10*8+32],ymm1
                vmovapd [rdx+r10*8+64],ymm2
                vmovapd [rdx+r10*8+96],ymm3
                vmovapd [rdx+r10*8+128],ymm4
                vmovapd [rdx+r10*8+160],ymm5
                vmovapd [rdx+r10*8+192],ymm6
                vmovapd [rdx+r10*8+224],ymm7
                vmovapd [rdx+r10*8+256],ymm8
                vmovapd [rdx+r10*8+288],ymm9
                vmovapd [rdx+r10*8+320],ymm10
                vmovapd [rdx+r10*8+352],ymm11
                vmovapd [rdx+r10*8+384],ymm12
                vmovapd [rdx+r10*8+416],ymm13
                vmovapd [rdx+r10*8+448],ymm14
                vmovapd [rdx+r10*8+480],ymm15
                add r10, 64
                cmp r10, rbx
                jl cicloV
                
                xor r10, r10
                xor r11, r11
                mov r12d, [rsi]
          ciclo:
                vxorpd ymm0, ymm0
                cmp r12, 0
                je fine
                
         ciclor:
                mov r13d, [rdi+r11*4]
                sal r13d, 3
                add r13, rdx
                vmovsd xmm1, [r13]
                vaddsd xmm0, xmm1
   
                inc r11
                dec r12
                cmp r12, 0
                jg ciclor
                
           fine:
                vmovsd [r9+r10*8], xmm0
                mov r12d, [rsi+r10*4+4]
                sub r12d, [rsi+r10*4]
                inc r10  
                cmp r10, rcx
                jl ciclo
                

		popaq						; ripristina i registri generali
		mov		rsp, rbp			; ripristina lo Stack Pointer
		pop		rbp				; ripristina il Base Pointer
		ret						; torna alla funzione C chiamante
