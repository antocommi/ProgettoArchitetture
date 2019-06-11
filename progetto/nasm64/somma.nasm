%include "sseutils64.nasm"

section .data
section .bss
section .text

global somma

; dim   	rdx
; dest  	rsi
; source	rdi

somma:
		push rbp
		mov	rbp, rsp
		;pushaq

        xor rcx, rcx
        sub rdx, 4
cicloQ: cmp rcx, rdx
        jg endQ
        vmovaps xmm0, [rsi+4*rcx]
        vmovaps xmm1, [rdi+4*rcx]
        vaddps xmm0, xmm1
        vmovaps [rsi+4*rcx], xmm0
        add rcx, 4
        jmp cicloQ
endQ:   add rdx, 4

cicloR: cmp rcx, rdx
        jge endR
        vmovss xmm0, [rsi+4*rcx]
        vmovss xmm1, [rdi+4*rcx]
        vaddss xmm0, xmm1
        vmovss [rsi+4*rcx], xmm0
        inc rcx
        jmp cicloR
endR:
       ;popaq
		mov	rsp, rbp
		pop	rbp
		ret


