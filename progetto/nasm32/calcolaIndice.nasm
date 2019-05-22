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
		push edi	;inizio
		mov eax, [ebp+i]	;i
		mov ebx, [ebp+j]	;j
		mov ecx, eax		;i
		dec eax				;i-1
		imul eax, ecx		;i*(i-1)
		shr eax, 1			;i*(i-1)/2
		add eax, ebx		;i*(i-1)/2+j
		pop	edi		;fine
		pop	esi
		pop	ebx
		mov	esp, ebp
		pop	ebp
		ret