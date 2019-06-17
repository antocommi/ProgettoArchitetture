%include "sseutils64.nasm"

extern distanza

section .data
max dd 1.79E+308
; f db '%f ', 10, 0
; d1 db '%ld ';, 10, 0
; d2 db 'i%ld j%ld ', 10, 0
; d3 db '%ld %ld %ld ', 10, 0
; d4 db '%ld %ld %ld %ld ', 10, 0
; break db 'breakpoint', 10, 0

section .bss

section .text

global calcolaPQ

; end equ 20		rcx
; startt equ 16		rdx
; partition equ 12	rsi
; data equ 8		rdi

source equ 0
dim_source equ 8
index equ 16
dest equ 24
index_rows equ 32
index_columns equ 36
n_centroidi equ 40
d equ 44

; punto1 equ -48		rdi
; punto2 equ -44		rsi
; dimensione equ -40	rdx

calcolaPQ:
		push rbp
		mov	rbp, rsp
		push rbx
		push r12
		;pushaq
		;pushad		;inizio

		mov r11, rdx;start
		mov r8, rsi	;partition
		mov r9, rdi ;data

		mov ebx, [r9+d]
		mov r12, rbx				;d
		sal r12, 2
		mov ebx, [r9+index_columns]	;m
		sal rbx, 2

		sub rcx, rdx
		mov rdx, rcx ;primo parametro

		mov rdi, r11
		sal rdi, 2
		add rdi, [r9+source];ind1=source+start

		mov r10, r8
		sal r10, 2
		add r10, [r9+index];ind=index+partition

		mov ecx, [r9+n_centroidi]
		mov r8, rcx

		vmovss xmm6, [max]

		mov ecx, [r9+dim_source]
forI:	cmp rcx, 0
		jle endI
		mov rsi, r11
		sal rsi, 2
		add rsi, [r9+dest];ind2=dest+start

; 		pushaq
;  		cmp r11, 0
;  		je boh
;  		mov rsi, rcx
;  		mov rdi, d1
;  		call printf
; boh:	popaq
		
		vmovss xmm5, xmm6	;max
		xor rax, rax
forJ:	cmp rax, r8
		jge endJ
		
;  		pushaq
;  		cmp r11, 0
;  		je boh
;  		mov rcx, rdx
;  		mov rdx, rsi
;  		mov rsi, rdi
;  		mov rdi, d3
;  		call printf
;  boh:	popaq

		call distanza

; 		pushaq
; 		cmp r11, 0
; 		je boh2
; 		mov rdx, rax
; 		mov rsi, rcx
; 		mov rdi, d2
; 		call printf
; boh2:	popaq

		vcomiss xmm0, xmm5
		jae endif

		vmovss xmm5, xmm0
		mov [r10], eax
endif:	add rsi, r12
		inc rax

		; pushaq
		; mov edi, break
		; call printf
		; popaq

		jmp forJ
endJ:	add r10, rbx
		add rdi, r12
		dec rcx

; 		pushaq
;  		cmp r11, 0
;  		je boh1
;  		mov rsi, rcx
;  		mov rdi, d1
;  		call printf
; boh1:	popaq

		jmp forI
endI:	;popad		;fine
		;popaq
		pop r12
		pop rbx
		mov	rsp, rbp
		pop	rbp
		ret