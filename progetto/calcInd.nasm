%include "sseutils.nasm"

section .data

section .bss

section .text

global calcolaIndice

j equ 8
i equ 12

calcolaIndice:
		push ebp
		mov	ebp, esp
		push ebx
		push esi
		push edi

		mov eax, [ebp+i]
		mov ebx, [ebp+j]
		mov ecx, eax
		dec eax
		imul eax, ecx
		shr eax, 1
		add eax, ebx

		pop	edi
		pop	esi
		pop	ebx
		mov	esp, ebp
		pop	ebp
		ret
