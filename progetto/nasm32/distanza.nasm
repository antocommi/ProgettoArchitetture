section .data
section .bss
section .text

global distanza

r equ 20
dimensione equ 16
punto2 equ 12
punto1 equ 8

distanza:
		push ebp
		mov	ebp, esp
		pushad	;inizio
		
		mov ecx, [ebp+punto1]
		mov edx, [ebp+punto2]
		mov ebx, [ebp+dimensione]

		xorps xmm0, xmm0
		xor esi, esi				;i=0

		sub ebx, 8
cicloQ:	cmp esi, ebx				;i < end-start
		jg somme
		
		movups xmm1, [ecx+4*esi]
		movups xmm3, [ecx+4*esi+16]
		movups xmm2, [edx+4*esi]
		movups xmm4, [edx+4*esi+16]

		subps xmm1, xmm2
		subps xmm3, xmm4
		mulps xmm1, xmm1
		mulps xmm3, xmm3
		addps xmm0, xmm1
		addps xmm0, xmm3
		
		add esi, 8
		jmp cicloQ
somme:	haddps xmm0, xmm0
		haddps xmm0, xmm0
		add ebx, 8

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
		popad		;fine
		mov	esp, ebp
		pop	ebp
		ret