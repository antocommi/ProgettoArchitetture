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
		push ebx
		push esi
		push edi	;inizio
		mov ecx, [ebp+punto1]
		mov edx, [ebp+punto2]
		mov ebx, [ebp+dimensione]
		sub ebx, 4
		xorps xmm0, xmm0
		xor esi, esi				;i=0
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