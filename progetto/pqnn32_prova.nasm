%include "sseutils.nasm"

section .data
section .bss
section .text

global calcolaIndice
global dist_eI

;calcolaIndice
j equ 8
i equ 12

;dist_eI
end equ 8
start equ 12
punto2 equ 16
punto1 equ 20
set equ 24
input equ 28

calcolaIndice:
		push ebp
		mov	ebp, esp
		push ebx
		push esi
		push edi	;inizio
		mov eax, [ebp+i]
		mov ebx, [ebp+j]
		mov ecx, eax
		dec eax
		imul eax, ecx
		shr eax, 1
		add eax, ebx
		pop	edi		;fine
		pop	esi
		pop	ebx
		mov	esp, ebp
		pop	ebp
		ret

dist_eI:
		push ebp
		mov	ebp, esp
		push ebx
		push esi
		push edi	;inizio
		xor esi, esi			;i=0
		mov ebx, [ebp+end]		;end
		sub ebx, [ebp+start]	;end-start
		sub ebx, 4				;end ciclo quoziente
		mov ecx, [ebp+input+28]	;ind=input->d
		mov edx, ecx			;ind2=input->d
		imul ecx, [ebp+punto1]	;ind=input->d*punto1
		add ecx, [ebp+set]		;ind=input->d*punto1+set
		add ecx, ebx			;ind=input->d*punto1+set+start
		imul edx, [ebp+punto2]	;ind2=input->d*punto2
		add edx, [ebp+input+8]	;ind2=input->d*punto2+input->ds
		add edx, ebx			;ind2=input->d*punto2+input->ds+start
		xor xmm1, xmm1
cicloQ:	cmp esi, ebx			;i < end-start
		jg somme
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
		jg endloop
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