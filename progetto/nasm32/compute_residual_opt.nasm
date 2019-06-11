%include "sseutils.nasm"

extern printf

section .data
break db 'breakpoint', 10, 0
section .bss
section .text

global compute_residual_opt

src equ 24
y equ 20
qc_i equ 16
res equ 12
input equ 8


qc equ 112
d equ 16


compute_residual_opt:
	push ebp
	mov	ebp, esp
	pushad    ;inizio

        mov eax, [ebp+input]
        mov ebx, [eax+d]
        mov ecx, ebx 
        imul ebx, [ebp+y]
        imul ebx, 4             ; poi modificare  on shift a sinistra
        add ebx, [ebp+src]

        imul ecx, [ebp+qc_i]
        imul ecx, 4
        add ecx, [eax+qc]

        mov edx, [ebp+res]

        xor esi, esi
        mov eax, [eax+d] ;eax=d
        sub eax, 16     
cicloQ: cmp esi, eax
        jg fineQ
        movaps xmm0, [ebx+4*esi]
        movaps xmm1, [ebx+4*esi+16]
        movaps xmm2, [ebx+4*esi+32]
        movaps xmm3, [ebx+4*esi+48]
        subps xmm0, [ecx+4*esi]
        subps xmm1, [ecx+4*esi+16]
        subps xmm2, [ecx+4*esi+32]
        subps xmm3, [ecx+4*esi+48]
        movaps [edx+4*esi], xmm0
        movaps [edx+4*esi+16], xmm1
        movaps [edx+4*esi+32], xmm2
        movaps [edx+4*esi+48], xmm3
        add esi,16
        jmp cicloQ
       
fineQ:  add eax, 16
cicloR1:cmp esi, eax
        jg fineR1
        movaps xmm0, [ebx+4*esi]
        subps xmm0, [ecx+4*esi]
        movaps [edx+4*esi], xmm0
        add esi, 4
        jmp cicloR1
fineR1: add eax, 4       
cicloR2:cmp esi, eax
        jg fine
        movss xmm0, [ebx+4*esi]
        subss xmm0, [ecx+4*esi]
        movss [edx+4*esi], xmm0 
        inc esi
        jmp cicloR2

       

fine:   popad ;fine
	mov esp, ebp
	pop ebp
	ret