%include "sseutils.nasm"

section .data
section .bss
section .text

global compute_residual_opt

y equ 20
qc_i equ 16
res equ 12
input equ 8

dataset equ 4
qc equ 112
d equ 16


compute_residual_opt:
		push ebp
		mov	ebp, esp
		push ebx
		push esi
		push edi    ;inizio

        mov eax, [ebp+input]
        mov ebx, [eax+d]
        mov ecx, ebx 
        imul ebx, [ebp+y]
        imul ebx, 4             ; poi modificare  on shift a sinistra
        add ebx, [eax+dataset]

        imul ecx, [ebp+qc_i]
        imul ecx, 4
        add ecx, [eax+qc]

        mov edx, [ebp+res]

        xor esi, esi
        mov eax, [eax+d] ;eax=d

ciclo:  movaps xmm0, [ebx+4*esi]
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
        cmp esi, eax
        jl ciclo

        pop	edi		;fine
		pop	esi
		pop	ebx
		mov	esp, ebp
		pop	ebp
		ret