%include "sseutils64.nasm"

extern distanza

section .data
d1 db '%d ', 10, 0
d5 db '%ld %ld %ld %ld %ld ', 10, 0
d3 db '%ld %ld %ld', 10, 0
f db '%f ', 10, 0
;break db 'breakpoint', 10, 0

section .bss

section .text

global calcolaFob

; end equ 24	r8
; start equ 20	rcx
; ipart equ 16	rdx		r11
; data equ 12	rsi		r10
; input equ 8	rdi		r9

;params
m equ 40

;kmeans_data
source equ 0
dim_source equ 8
index equ 16
dest equ 24
d equ 44

; punto1 equ -48		rdi
; punto2 equ -44		rsi
; dimensione equ -40	rdx


calcolaFob:
		push rbp
		mov	rbp, rsp
		;push rax
		push rbx
		push r12
		push r13
		;pushaq
		;pushad		;inizio

		mov r9, rdi		;input
		mov r10, rsi	;data
		mov r11, rdx	;ipart

		mov rdx, r8
		sub rdx, rcx	;parametro dimensione

		mov rsi, rcx	;ind2=start
		sal rsi, 2
		add rsi, [r10+source];ind2=source+start

		mov r12, rcx	;ind=start
		sal r12, 2
		add r12, [r10+dest];ind=dest+start

		mov r13, r11	;ipart
		sal r13, 2
		add r13, [r10+index];index+ipart
		
		mov eax, [r10+d]
		sal eax, 2
		mov ecx, [r9+m]
		sal ecx, 2

		vxorps xmm10, xmm10
		mov ebx, [r10+dim_source]

		sub rbx, 4
cicloQ:	cmp rbx, 0
		jl somme

		mov edi, [r13]
		imul rdi, rax
		add rdi, r12
		call distanza
		add rsi, rax
		add r13, rcx
		vmovss xmm11, xmm0

		mov edi, [r13]
		imul rdi, rax
		add rdi, r12
		call distanza
		add rsi, rax
		add r13, rcx
		vinsertps xmm11, xmm0, 00010000b

		mov edi, [r13]
		imul rdi, rax
		add rdi, r12
		call distanza
		add rsi, rax
		add r13, rcx
		vinsertps xmm11, xmm0, 00100000b

		mov edi, [r13]
		imul rdi, rax
		add rdi, r12
		call distanza
		add rsi, rax
		add r13, rcx
		vinsertps xmm11, xmm0, 00110000b

		;vmulps xmm11, xmm11
		vaddps xmm10, xmm11

		sub rbx, 4
		jmp cicloQ
somme:	add rbx, 4
cicloR:	cmp rbx, 0
		jle fine

		mov edi, [r13]
		imul rdi, rax
		add rdi, r12
		call distanza
		add rsi, rax
		add r13, rcx

		;vmulss xmm0, xmm0
		vaddss xmm10, xmm0

		dec rbx
		jmp cicloR
fine:	vmovss xmm0, xmm10
		;popaq
		;popad		;fine
		pop r13
		pop r12
		pop rbx
		;pop rax
		mov	rsp, rbp
		pop	rbp
		ret