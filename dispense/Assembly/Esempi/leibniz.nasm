; --------------------------------------------------
; Esempi di utilizzo delle istruzioni SSE
; --------------------------------------------------
; Formula di Leibniz per il calcolo di pi greco
; --------------------------------------------------
; F. Angiulli
;

;
; Per eseguire:
;
;     ./run leibniz
;
; Per visualizzare un maggior numero di cifre 
; significative nella stampa a video occorre 
; modificare la linea 22 del file sseutils.nasm 
; ad esempio come segue
;
; dmask:	db		'%18.16f ',0
;

%include "sseutils.nasm"

section .data			; Sezione contenente dati inizializzati

align 16 
inizio:		dq		3.0, 1.0
align 16
quattro:	dq		4.0, 4.0
align 16
uno:		dq		1.0, 1.0
pi:			dq		3.14159265358979

section .bss			; Sezione contenente dati non inizializzati

ris:		resq	1

section .text			; Sezione contenente il codice macchina

global	main

main:	start

		xorpd		xmm0, xmm0
		movapd		xmm1, [inizio]
		movapd		xmm2, [quattro]
		movapd		xmm3, [uno]
		mov			ecx, 5000000
ciclo:	
		movapd		xmm4, xmm3
		divpd		xmm4, xmm1
		addsubpd	xmm0, xmm4
		addpd		xmm1, xmm2
		dec		ecx
		jnz		ciclo
		haddpd		xmm0, xmm0
		addsd		xmm0, xmm0
		addsd		xmm0, xmm0
		movsd		[ris], xmm0

		printsd		ris
		printsd		pi

		stop
