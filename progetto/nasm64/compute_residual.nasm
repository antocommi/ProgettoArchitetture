%include "sseutils64.nasm"

section .data
section .bss
section .text

global compute_residual

; src equ 24	r8
; y equ 20    rcx
; qc_i equ 16 rdx
; res equ 12  rsi
; input equ 8 rdi


; qc equ 112
; d equ 16

qc equ 160
d equ 28


compute_residual:
		push rbp
		mov	rbp, rsp
		;pushaq
		
		mov eax, [rdi+d]
		mov r9, rax
        mov r10, r9 
		imul r9, rcx
		imul r9, 4             ; poi modificare  on shift a sinistra
        add r9, r8

		imul r10, rdx
		imul r10,4
		add r10, [rdi+qc]
       
        xor r11, r11
        mov edi, eax ;rdi=d
        sub rdi, 64     ;????
		
cicloQ1: cmp r11, rdi
        jg fineQ1
		vmovups ymm0, [r9+r11*4]
        vmovups ymm1, [r9+r11*4+32]
        vmovups ymm2, [r9+r11*4+64]
        vmovups ymm3, [r9+r11*4+96]        
        vmovups ymm4, [r9+r11*4+128]
        vmovups ymm5, [r9+r11*4+160]
        vmovups ymm6, [r9+r11*4+192]
        vmovups ymm7, [r9+r11*4+224]

        vsubps ymm0, [r10+r11*4]
        vsubps ymm1, [r10+r11*4+32] 
        vsubps ymm2, [r10+r11*4+64] 
        vsubps ymm3, [r10+r11*4+96] 
        vsubps ymm4, [r10+r11*4+128]         
        vsubps ymm5, [r10+r11*4+160] 
        vsubps ymm6, [r10+r11*4+192] 
        vsubps ymm7, [r10+r11*4+224]
		
        vmovups [rsi+r11*4],ymm0 
        vmovups [rsi+r11*4+32],ymm1
        vmovups [rsi+r11*4+64],ymm2 
        vmovups [rsi+r11*4+96],ymm3 
        vmovups [rsi+r11*4+128],ymm4 
        vmovups [rsi+r11*4+160],ymm5 
        vmovups [rsi+r11*4+192],ymm6 

		add r11,64 ;????
        jmp cicloQ1
fineQ1:  add rdi, 56
cicloR1:cmp r11, rdi
        jg fineR1
		vmovups ymm0, [r9+r11*4]
        vsubps ymm0, [r10+r11*4]
        vmovups [rsi+r11*4],ymm0 
        add r11,8 ;????
        jmp cicloQ1
fineR1:  add rdi, 8
cicloR2: cmp r11, rdi
        jge fine
        vmovss xmm0, [r9+4*r11]
        vsubss xmm0, [r10+4*r11]
        vmovss [rsi+4*r11], xmm0 
        inc r11
        jmp cicloR2

fine:	;popaq						; ripristina i registri generali
		mov		rsp, rbp			; ripristina lo Stack Pointer
		pop		rbp				; ripristina il Base Pointer
		ret						; torna alla funzione C chiamante