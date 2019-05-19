%include "sseutils.nasm"

section .data
section .bss
section .text

global dist_eI

end equ 8
startt equ 12
punto2 equ 16
punto1 equ 20
set equ 24
input equ 28

in_d equ 28
in_ds equ 8

dist_eI:
		push ebp
		mov	ebp, esp
		push ebx
		push esi
		push edi	;inizio
		xor esi, esi				;i=0
		mov ebx, [ebp+end]			;end
		sub ebx, [ebp+startt]		;end-start
		sub ebx, 4					;end ciclo quoziente
		mov ecx, [ebp+input+in_d]	;ind=input->d
		mov edx, ecx				;ind2=input->d
		imul ecx, [ebp+punto1]		;ind=input->d*punto1
		add ecx, [ebp+set]			;ind=input->d*punto1+set
		add ecx, ebx				;ind=input->d*punto1+set+start
		imul edx, [ebp+punto2]		;ind2=input->d*punto2
		mov edi, [ebp+input+in_ds]	;input->ds
		add edx, [edi]				;ind2=input->d*punto2+input->ds
		add edx, ebx				;ind2=input->d*punto2+input->ds+start
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
endloop:pop	edi		;fine
		pop	esi
		pop	ebx
		mov	esp, ebp
		pop	ebp
		ret