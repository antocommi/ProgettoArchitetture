extern calcolaIndice

section .data
section .bss
section .text

global dist_matrix

ipart equ 20
centroide2 equ 16
centroide1 equ 12
input equ 8

m equ 28
distanze_simmetriche equ 96
zero equ 124

dist_matrix:
		push ebp
		mov	ebp, esp
		push ebx
		push ecx
		push edx
		push esi
		push edi
		;pushad	;inizio
		mov esi, [ebp+centroide1]
		mov edi, [ebp+centroide2]
		cmp esi, edi
		jne else
		mov ebx, [ebp+input]
		mov eax, [ebx+zero]
		jmp fine
else:	;cmp esi, edi
		jg else2
		push esi
		push edi
		jmp cal
else2:	push edi
		push esi
cal:	call calcolaIndice
		add esp, 8
		mov ebx, [ebp+input]
		imul eax, [ebx+m]
		add eax, [ebp+ipart]
		imul eax, 4
		add eax, [ebx+distanze_simmetriche]
fine:	;popad		;fine
		pop edi
		pop esi
		pop edx
		pop ecx
		pop ebx
		mov	esp, ebp
		pop	ebp
		ret