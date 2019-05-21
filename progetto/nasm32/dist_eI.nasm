%include "sseutils.nasm"

section .data
section .bss
section .text

global dist_eI

r equ 32
end equ 28
startt equ 24
punto2 equ 20
punto1 equ 16
set equ 12
input equ 8

in_d equ 16
in_ds equ 4

dist_eI:
		push ebp
		mov	ebp, esp
		push ebx
		push esi
		push edi	;inizio
		mov ebx, [ebp+end]			;end
		mov edi, [ebp+input]		;input
		mov ecx, [edi+in_d]			;ind=input->d
		mov edx, ecx				;ind2=input->d
		imul ecx, [ebp+punto1]		;ind=input->d*punto1
		add ecx, [ebp+set]			;ind=input->d*punto1+set
		add ecx, [ebp+startt]		;ind=input->d*punto1+set+start
		imul edx, [ebp+punto2]		;ind2=input->d*punto2
		add edx, [edi+in_ds]		;ind2=input->d*punto2+input->ds
		add edx, [ebp+startt]		;ind2=input->d*punto2+input->ds+start
		sub ebx, [ebp+startt]		;end-start
		sub ebx, 4					;end ciclo quoziente
		xorps xmm0, xmm0
		xor esi, esi				;i=0
cicloQ:	cmp esi, ebx				;i < end-start
		jge somme
		movaps xmm1, [ecx+4*esi]
		subps xmm1, [edx+4*esi]
		mulps xmm1, xmm1
		addps xmm0, xmm1
		add esi, 4
		jmp cicloQ
somme:	haddps xmm0, xmm0
		haddps xmm0, xmm0
		add ebx, 4
cicloR:	cmp esi, ebx
		jge endloop
		movss xmm1, [ecx+4*esi]
		subss xmm1, [edx+4*esi]
		mulss xmm1, xmm1
		addss xmm0, xmm1
		inc esi
		jmp cicloR
endloop:mov eax, [ebp+r];
		movss [eax], xmm0
		pop	edi		;fine
		pop	esi
		pop	ebx
		mov	esp, ebp
		pop	ebp
		ret