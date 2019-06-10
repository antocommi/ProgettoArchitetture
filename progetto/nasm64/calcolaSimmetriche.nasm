%include "sseutils64.nasm"

extern dist

section .data
d1 db '%ld ', 10, 0
d5 db '%ld %ld %ld %ld %ld ', 10, 0
d3 db '%ld %ld %ld', 10, 0
f db '%f ', 10, 0
break db 'breakpoint', 10, 0

section .bss

section .text

global calcolaSimmetriche

; query	rdx
; ind	rsi
; input	rdi

n equ 24
query_pq equ 112

; punto2 equ 20		rcx
; punto1 equ 16		rdx
; quantizer equ 12	rsi
; input equ 8		rdi

calcolaSimmetriche:
		push rbp
		mov	rbp, rsp
		;pushaq

		mov r8, rsi;ind
		mov rsi, [rdi+query_pq]

		mov eax, [rdi+n]
		xor rcx, rcx

		sub rax, 8

cicloQ:	cmp rcx, rax
		jg endQ

		call dist
		vmovss xmm11, xmm0
		inc rcx
		call dist
		vinsertps xmm11, xmm0, 00010000b
		inc rcx
		call dist
		vinsertps xmm11, xmm0, 00100000b
		inc rcx
		call dist
		vinsertps xmm11, xmm0, 00110000b
		inc rcx

		;sposta xmm11 nella parte alta di ymm11
		;vperm2f128 ymm11, ymm11, ymm11, 00000000b

		call dist
		vmovss xmm12, xmm0
		inc rcx
		call dist
		vinsertps xmm12, xmm0, 00010000b
		inc rcx
		call dist
		vinsertps xmm12, xmm0, 00100000b
		inc rcx
		call dist
		vinsertps xmm12, xmm0, 00110000b
		inc rcx

		vperm2f128 ymm13, ymm11, ymm12, 00100000b
		;vperm2f128 ymm11, ymm11, ymm13, 00100000b

		vmovaps [r8+4*rcx-32], ymm13

		jmp cicloQ
endQ:	add rax, 4

cicloQ2:cmp rcx, rax
		jg endQ2

		call dist
		vmovss xmm11, xmm0
		inc rcx
		call dist
		vinsertps xmm11, xmm0, 00010000b
		inc rcx
		call dist
		vinsertps xmm11, xmm0, 00100000b
		inc rcx
		call dist
		vinsertps xmm11, xmm0, 00110000b
		inc rcx

		vmovaps [r8+4*rcx-16], xmm11

		jmp cicloQ2
endQ2:	add rax, 4

cicloR:	cmp rcx, rax
		jge endR

		call dist
		vmovss [r8+4*rcx], xmm0
		inc rcx

		jmp cicloR
endR:   ;popaq
		mov	rsp, rbp
		pop	rbp
		ret








