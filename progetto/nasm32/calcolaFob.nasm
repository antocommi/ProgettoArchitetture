extern distanza

section .data
section .bss
temp: resd 1

section .text

global calcolaFob


r equ 28
end equ 24
start equ 20
ipart equ 16
data equ 12
input equ 8

;kmeans_data
source equ 0
dim_source equ 4
index equ 8
dest equ 12
index_columns equ 20
d equ 28

punto1 equ -48
punto2 equ -44
dimensione equ -40
cr equ -36


calcolaFob:
		push ebp
		mov	ebp, esp
		pushad		;inizio
		sub esp, 16

		mov dword [ebp+cr], temp			;primo parametro distanza
		mov eax, [ebp+end]
		mov ebx, [ebp+start]
		sub eax, ebx
		mov [ebp+dimensione], eax			;secondo parametro distanza
		mov ecx, [ebp+data]					;ecx=data
		mov edx, [ebp+input]				;edx=input
		imul ebx, 4							;solito problema
		
		mov eax, ebx
		add ebx, [ecx+source]
		mov [ebp+punto2], ebx				;terzo parametro distanza
		mov esi, [ecx+dim_source]			;ciclo

		mov ebx, [ebp+ipart]
		imul ebx, 4
		add ebx, [ecx+index]				;ebx=index[ipart]
		add eax, [ecx+dest]					;eax=dest+start

		mov edx, [ecx+index_columns]		;edx=m
		imul edx, 4
		mov ecx, [ecx+d]					;ecx=d
		imul ecx, 4
		xorps xmm7, xmm7
		
		mov edi, [ebx]
		imul edi, ecx
		add edi, eax
		mov [ebp+punto1], edi

ciclo:	cmp esi, 0
		jle fine

		call distanza
		
		addss xmm7, xmm0
		add [ebp+punto2], edx		;aggiornamento terzo parametro

		add ebx, edx
		mov edi, [ebx]
		imul edi, ecx
		add edi, eax
		mov [ebp+punto1], edi		;aggiornamento quarto parametro

		dec esi
		jmp ciclo
fine:	mov eax, [ebp+r]
		movss [eax], xmm7
		add esp, 16
		popad		;fine
		mov	esp, ebp
		pop	ebp
		ret