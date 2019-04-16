; --------------------------------------------------
; Esempio di utilizzo delle istruzioni SSE
;
; Dati due vettori a ed x di 4 single precision, 
; calcola il vettore b i cui elementi sono dati da:
;
;     b[i] = a[i]*x[i]^2 + 2*x[i]/a[i] - 3
;
; --------------------------------------------------
; F. Angiulli
;

;
; Per eseguire:
;
;     ./run esempio1
;

%include "sseutils.nasm"

section .data			; Sezione contenente dati inizializzati

amsg:	db		'a = ',0
xmsg:	db		'x = ',0
bmsg:	db		'b = ',0

align 16 
a:		dd		1.5, 2.3, 1.1, -3.5
align 16
x:		dd		1.0, -3.0, 5.0, 4.0
align 16
tre:	dd		3.0, 3.0, 3.0, 3.0

section .bss			; Sezione contenente dati non inizializzati

alignb 16
b:		resd	4

section .text			; Sezione contenente il codice macchina

global	main

main:	start

		movaps		xmm0, [x]
		movaps		xmm1, xmm0
		mulps		xmm0, xmm0
		movaps		xmm2, [a]
		mulps		xmm0, xmm2
		addps		xmm1, xmm1
		divps		xmm1, xmm2
		addps		xmm0, xmm1
		subps		xmm0, [tre]
		movaps		[b], xmm0

		prints		amsg
		printps		a, 1
		prints		xmsg
		printps		x, 1
		prints		bmsg
		printps		b, 1

		stop
