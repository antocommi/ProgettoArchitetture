extern dist_matrix
extern printf

section .data
;d1 db '%d ', 10, 0
f db '%f ', 10, 0
;break db 'breakpoint', 10, 0
section .bss
section .text

global dist

r equ 24
punto2 equ 20
punto1 equ 16
quantizer equ 12
input equ 8

m equ 28
pq equ 84
distanze_simmetriche equ 96

inputc equ -48
c1 equ -44
c2 equ -40
i equ -36

dist:
		push ebp
		mov	ebp, esp
		pushad		;inizio
		sub esp, 16

		mov ebx, [ebp+input]		;input
		mov ecx, [ebx+m]			;c1=m
		mov edx, ecx				;c2=m
		mov edi, ecx
		imul ecx, [ebp+punto1]		;c1=m*punto1
		imul edx, [ebp+punto2]		;c2=m*punto2
		imul ecx, 4					;solito problema
		imul edx, 4					;solito problema
		add ecx, [ebp+quantizer]	;c1=m*punto1+quantizer
		add edx, [ebx+pq]			;c2=m*punto2+pq
		mov [ebp+inputc], ebx		;  primo parametro
		sub edi, 4

		xorps xmm0, xmm0
		xor esi, esi
		mov dword [ebp+i], 0
cicloQ:	cmp esi, edi
		jg endQ
		;-----------------------------------------
		mov eax, [ecx+4*esi]
		mov [ebp+c1], eax			;  secondo parametro		
		mov eax, [edx+4*esi]
		mov [ebp+c2], eax			;  terzo parametro
		call dist_matrix
		pinsrd xmm1, [eax], 0
		inc dword [ebp+i]			;  quarto parametro

		mov eax, [ecx+4*esi+4]
		mov [ebp+c1], eax			;  secondo parametro		
		mov eax, [edx+4*esi+4]
		mov [ebp+c2], eax			;  terzo parametro
		call dist_matrix
		pinsrd xmm1, [eax], 1
		inc dword [ebp+i]			;  quarto parametro

		mov eax, [ecx+4*esi+8]
		mov [ebp+c1], eax			;  secondo parametro		
		mov eax, [edx+4*esi+8]
		mov [ebp+c2], eax			;  terzo parametro
		call dist_matrix
		pinsrd xmm1, [eax], 2
		inc dword [ebp+i]			;  quarto parametro

		mov eax, [ecx+4*esi+12]
		mov [ebp+c1], eax			;  secondo parametro		
		mov eax, [edx+4*esi+12]
		mov [ebp+c2], eax			;  terzo parametro
		call dist_matrix
		pinsrd xmm1, [eax], 3
		inc dword [ebp+i]			;  quarto parametro

		addps xmm0, xmm1
		;-----------------------------------------
		add esi, 4
		jmp cicloQ

endQ:	haddps xmm0, xmm0
		haddps xmm0, xmm0
		add edi, 4

cicloR:	cmp esi, edi
		jge fine
		mov eax, [ecx+4*esi]
		mov [ebp+c1], eax			;  secondo parametro		
		mov eax, [edx+4*esi]
		mov [ebp+c2], eax			;  terzo parametro
		mov [ebp+i], esi			;  quarto parametro
		call dist_matrix
		movss xmm1, [eax]
		addss xmm0, xmm1
		inc esi
		jmp cicloR
fine:	mov eax, [ebp+r]
		movss [eax], xmm0
		add esp, 16
		popad		;fine
		mov	esp, ebp
		pop	ebp
		ret