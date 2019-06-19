%include "sseutils64.nasm"

extern dist_matrix
;extern printf

section .data
d1 db '%ld ', 10, 0
d4 db '%ld %ld %ld %ld', 10, 0
f db '%lf ', 10, 0
;break db 'breakpoint', 10, 0
section .bss
section .text

global dist

; punto2 equ 20		rcx
; punto1 equ 16		rdx
; quantizer equ 12	rsi
; input equ 8		rdi

m equ 40
pq equ 104
distanze_simmetriche equ 128

; ipart equ 20		rcx
; centroide2 equ 16	rdx
; centroide1 equ 12	rsi
; input equ 8		rdi

dist:
		push rbp
		mov	rbp, rsp
		push rsi
		push rdx
		push rcx
		push rax
		push rbx
		;pushaq
		;push rcx
		;pushad		;inizio

		mov r9, rdx		;punto1
		mov r10, rcx	;punto2
		mov ebx , [rdi+m]
		imul r9, rbx
		imul r10, rbx
		sal r9, 2
		sal r10, 2
		add r9, rsi
		add r10, [rdi+pq]

		vxorps xmm0, xmm0
		xor rcx, rcx
		sub rbx, 4
cicloQ:	cmp rcx, rbx
		jg fineQ

		mov esi, [r9+4*rcx]
		mov edx, [r10+4*rcx]
		call dist_matrix
		vpinsrd xmm1, [rax], 0
		inc rcx

		mov esi, [r9+4*rcx]
		mov edx, [r10+4*rcx]
		call dist_matrix
		vpinsrd xmm1, [rax], 1
		inc rcx
		
		mov esi, [r9+4*rcx]
		mov edx, [r10+4*rcx]
		call dist_matrix
		vpinsrd xmm1, [rax], 2
		inc rcx

		mov esi, [r9+4*rcx]
		mov edx, [r10+4*rcx]
		call dist_matrix
		vpinsrd xmm1, [rax], 3
		inc rcx

		vaddps xmm0, xmm1

		jmp cicloQ
fineQ:	vhaddps xmm0, xmm0
		vhaddps xmm0, xmm0
		add rbx, 4
		
cicloR:	cmp rcx, rbx
		jg fine
		mov esi, [r9+4*rcx]
		mov edx, [r10+4*rcx]
		call dist_matrix
		
		vaddss xmm0, [rax]

		inc rcx
		jmp cicloR

fine:	;popad		;fine
		;pop rcx
		;popaq
		pop rbx
		pop rax
		pop rcx
		pop rdx
		pop rsi
		mov	rsp, rbp
		pop	rbp
		ret