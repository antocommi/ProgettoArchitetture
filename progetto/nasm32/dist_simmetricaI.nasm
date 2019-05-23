%include "sseutils.nasm"

section .data
section .bss
section .text

global dist_simmetricaI

r equ 28
end equ 24
startt equ 20
centroide2 equ 16
centroide1 equ 12
input equ 8

in_d equ 16
in_codebook equ 92

dist_simmetricaI:
		push ebp
		mov	ebp, esp
		push ebx
		push esi
		push edi	;inizio
		xor esi, esi				;i=0
		mov ebx, [ebp+end]			;end
		sub ebx, [ebp+startt]		;end-start
		sub ebx, 4					;end ciclo quoziente
		mov edi, [ebp+input]		;input
		mov ecx, [edi+in_d]			;ind=input->d
		mov edx, ecx				;ind2=input->d
		imul ecx, [ebp+centroide1]	;ind=input->d*centroide1
		add ecx, [edi+in_codebook]	;ind=input->d*centroide1+codebook
		add ecx, [ebp+startt]		;ind=input->d*centroide1+codebook+start
		imul edx, [ebp+centroide2]	;ind2=input->d*centroide2
		add edx, [edi+in_codebook]	;ind2=input->d*centroide2+codebook
		add edx, [ebp+startt]		;ind2=input->d*centroide2+codebook+start
		xorps xmm1, xmm1
cicloQ:	cmp esi, ebx				;i < end-start
		jge somme
		movaps xmm0, [ecx+4*esi]
		subps xmm0, [edx+4*esi]
		mulps xmm0, xmm0
		addps xmm1, xmm0
		add esi, 4
		jmp cicloQ
somme:	haddps xmm1, xmm1
		haddps xmm1, xmm1
		add ebx, 4
cicloR:	cmp esi, ebx
		jge endloop
		movss xmm0, [ecx+4*esi]
		subss xmm0, [edx+4*esi]
		mulss xmm0, xmm0
		addss xmm1, xmm0
		inc esi
		jmp cicloR
endloop:mov eax, [ebp+r]
		movss [eax], xmm0
		pop	edi		;fine
		pop	esi
		pop	ebx
		mov	esp, ebp
		pop	ebp
		ret