%include "sseutils.nasm"

section .data
section .bss
section .text

global compute_residual

y equ 20
qc_i equ 16
res equ 12
input equ 8

dataset equ 4
qc equ 112
d equ 16


compute_residual:
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
        mov eax, [eax+d]
        sub eax, 4
cicloQ: cmp esi, eax
        jg fineQ
        movaps xmm0, [ebx+4*esi]
        subps xmm0, [ecx+4*esi]
        movaps [edx+4*esi], xmm0
        add esi, 4
        jmp cicloQ
fineQ:  add eax, 4
cicloR: cmp esi, eax
        jg fineR
        movss xmm0, [ebx+4*esi]
        subss xmm0, [ecx+4*esi]
        movss [edx+4*esi], xmm0
        inc esi
        jmp cicloR

fineR:  pop	edi		;fine
		pop	esi
		pop	ebx
		mov	esp, ebp
		pop	ebp
		ret