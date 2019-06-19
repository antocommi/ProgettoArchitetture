%include "sseutils64.nasm"

extern calcolaIndice

section .data
d1 db '%d ', 10, 0
;f db '%f ', 10, 0
;break db 'breakpoint', 10, 0
section .bss
section .text

global dist_matrix

; ipart equ 20		rcx
; centroide2 equ 16	rdx
; centroide1 equ 12	rsi
; input equ 8		rdi

m equ 40
distanze_simmetriche equ 128
zero equ 184

dist_matrix:
		push rbp
		mov	rbp, rsp
		push rbx
		push rdi
		;pushaq
		;push rdx
		;pushad	;inizio

		mov rbx, rdi	;input
		cmp rsi, rdx
		jne else
		mov rax, [rbx+zero]
		jmp fine
else:	;cmp esi, edi
		jg else2
		mov rdi, rdx
		jmp cal
else2:	mov rdi, rsi
		mov rsi, rdx
cal:	call calcolaIndice

		xor rdx, rdx
		mov edx, [rbx+m]
		imul rax, rdx
		add rax, rcx
		sal rax, 2
		add rax, [rbx+distanze_simmetriche]

fine:	;popad		;fine
		;pop rdx
		;popaq
		pop rdi
		pop rbx
		mov	rsp, rbp
		pop	rbp
		ret