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

 %include "sseutils.nasm"

section .data			; Sezione contenente dati inizializzati

;uno:		db		'1.0',0
;
;align 16
;inizio:		dd		1.0, 2.0, 3.0, 4.0

section .bss			; Sezione contenente dati non inizializzati

;alignb 16
ri:		resd	1
r:              resq    1  
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
extern printf

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

global prod

P		equ		8
v               equ             12
dim             equ             16
n               equ             20
temp            equ             24

prod:
		push		ebp; salva il Base Pointer
		mov		ebp, esp; il Base Pointer punta al Record di Attivazione corrente
		push		ebx; salva i registri da preservare
		push		esi
		push		edi
		
		mov		eax, [ebp+P]
                mov             ebx, [ebp+v]

                xor esi, esi
                xor ecx, ecx
         cicloi:xor edi, edi
                mov edx, [ebp+dim]
                xorpd xmm4, xmm4
                xorpd xmm5, xmm5
                xorpd xmm6, xmm6
                xorpd xmm7, xmm7
         cicloj:
                movapd xmm0, [eax+ecx*8]
                movapd xmm1, [eax+ecx*8+16] 
                movapd xmm2, [eax+ecx*8+32] 
                movapd xmm3, [eax+ecx*8+48]
                mulpd xmm0, [ebx+edi*8]
                mulpd xmm1, [ebx+edi*8+16]
                mulpd xmm2, [ebx+edi*8+32]
                mulpd xmm3, [ebx+edi*8+48]
                addpd xmm4, xmm0
                addpd xmm5, xmm1
                addpd xmm6, xmm2
                addpd xmm7, xmm3
                add ecx, 8
                add edi,8
                cmp edi, edx
                jl cicloj
                addpd xmm4, xmm5
                addpd xmm4, xmm6
                addpd xmm4, xmm7
                haddpd xmm4, xmm4
                ;movhlps xmm1, xmm4
                ;addsd xmm4, xmm1
                mov eax, [ebp+temp]
                movsd [eax+esi*8], xmm4
                mov eax, [ebp+P]
                inc esi
                mov edx, [ebp+n]
                cmp esi, edx
                jl cicloi
                
        
      

		pop	edi; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp; ripristina lo Stack Pointer
		pop	ebp; ripristina il Base Pointer
		ret	   ; torna alla funzione C chiamante
