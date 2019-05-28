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
		mov ebx, [ebp+end]			;end
		sub ebx, [ebp+startt]		;end-start
		sub ebx, 4					;end ciclo quoziente
		mov edi, [ebp+input]		;input
		mov ecx, [edi+in_d]			;ind=input->d
		mov edx, ecx				;ind2=input->d
		imul ecx, [ebp+centroide1]	;ind=input->d*centroide1
		imul edx, [ebp+centroide2]	;ind2=input->d*centroide2
		add ecx, [ebp+startt]		;ind=input->d*centroide1+start
		add edx, [ebp+startt]		;ind2=input->d*centroide2+start

		imul ecx, 4					;gli indici usati in c vanno moltiplicati per 4 per farli diventare indirizzi,
		imul edx, 4					;in quanto ogni elemento occupa 4 celle
		
		add ecx, [edi+in_codebook]	;ind=input->d*centroide1+start+codebook
		add edx, [edi+in_codebook]	;ind2=input->d*centroide2+start+codebook
		xor esi, esi				;i=0
		xorps xmm0, xmm0
cicloQ:	cmp esi, ebx				;i < end-start
		jg somme
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
endloop:mov eax, [ebp+r]
		movss [eax], xmm0
		pop	edi		;fine
		pop	esi
		pop	ebx
		mov	esp, ebp
		pop	ebp
		ret