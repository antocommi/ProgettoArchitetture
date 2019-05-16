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
		mul eax, ecx
		div eax, 2
		add eax, ebx
		push eax
		
		pop	edi
		pop	esi
		pop	ebx
		mov	esp, ebp
		pop	ebp
		ret