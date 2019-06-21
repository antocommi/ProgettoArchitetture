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

        sub rdx, 8
cicloQ: cmp rcx, rdx
        jg endQ

        vmovups ymm0, [rsi+4*rcx]
        vmovups ymm1, [rdi+4*rcx]
        vaddps ymm0, ymm1
        vmovups [rsi+4*rcx], ymm0
        add rcx, 8
        jmp cicloQ
endQ:   add rdx, 8

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


