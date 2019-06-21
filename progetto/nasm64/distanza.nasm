section .data
section .bss
section .text

global distanza

; dimensione equ 16	rdx
; punto2 equ 12		rsi
; punto1 equ 8		rdi

distanza:
		push rbp
		mov	rbp, rsp
		push rbx
		;pushaq
		;pushad	;inizio

		vxorps ymm0, ymm0
		xor rbx, rbx

		sub rdx, 32
cicloQ:	cmp rbx, rdx				;i < end-start
		jg endQ

		vmovups ymm1, [rdi+4*rbx]
		vmovups ymm3, [rdi+4*rbx+32]
		vmovups ymm5, [rdi+4*rbx+64]
		vmovups ymm7, [rdi+4*rbx+96]
		vmovups ymm2, [rsi+4*rbx]
		vmovups ymm4, [rsi+4*rbx+32]
		vmovups ymm6, [rsi+4*rbx+64]
		vmovups ymm8, [rsi+4*rbx+96]

		vsubps ymm1, ymm2
		vsubps ymm3, ymm4
		vsubps ymm5, ymm6
		vsubps ymm7, ymm8

		vmulps ymm1, ymm1
		vmulps ymm3, ymm3
		vmulps ymm5, ymm5
		vmulps ymm7, ymm7

		vaddps ymm0, ymm1
		vaddps ymm0, ymm3
		vaddps ymm0, ymm5
		vaddps ymm0, ymm7
		
		add rbx, 32
		jmp cicloQ
endQ:	add rdx, 24
cicloQR:cmp rbx, rdx
		jg somme

		vmovups ymm1, [rdi+4*rbx]
		vmovups ymm2, [rsi+4*rbx]

		vsubps ymm1, ymm2
		vmulps ymm1, ymm1
		vaddps ymm0, ymm1

		add rbx, 8
		jmp cicloQR
somme:	add rdx, 8
		vhaddps ymm0, ymm0
		vhaddps ymm0, ymm0
		vperm2f128 ymm1, ymm0, ymm0, 1
		vaddps ymm0, ymm1
		
cicloR:	cmp rbx, rdx
		jge endloop
		vmovss xmm1, [rdi+4*rbx]
		vsubss xmm1, [rsi+4*rbx]
		vmulss xmm1, xmm1
		vaddss xmm0, xmm1
		inc rbx
		jmp cicloR
endloop:;sqrtss xmm0, xmm0
		;popad		;fine
		;popaq
		pop rbx
		mov	rsp, rbp
		pop	rbp
		ret