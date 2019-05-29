extern distanza

extern printf

section .data
max dd 1.79E+308
d1 db '%d '
;d2 db '%d %d '
;break db 'breakpoint', 10, 0

section .bss
temp: resd 1

section .text

global calcolaPQ

end equ 16
startt equ 12
data equ 8

source equ 0
dim_source equ 4
index equ 8
dest equ 12
index_rows equ 16
index_columns equ 20
n_centroidi equ 24
d equ 28

calcolaPQ:
		push ebp
		mov	ebp, esp
		push ebx
		push esi
		push edi	;inizio
		
		mov edx, [ebp+data]				;edx=data
		mov eax, [ebp+startt]			;ind=start
		imul eax, [edx+index_columns]	;ind=start*index_columns
		mov ecx, dword [edx+d]
		cdq
		idiv ecx				;ind=start*index_columns/d
		mov edx, [ebp+data]				;edx=data
		imul eax, 4						; solito problema  poi modificare con shift
		add eax, [edx+index] 			;ind=start*index_columns/d+index

		push temp						;  primo parametro funzione dist
		mov ecx, [ebp+end]
		sub ecx, [ebp+startt]
		push ecx						;  secondo parametro funzione dist
		mov ecx, [ebp+startt]			;ind1=start
		imul ecx, 4						; solito problema  poi modificare con shift
		add ecx, [edx+source]			;ind1=start+source
		push ecx						;  terzo parametro funzione dist
		mov ebx, [edx+d]
		imul ebx, 4

		mov esi, [edx+dim_source]		;end loop i
forI:	cmp esi, 0
		jle endI						;controllare condizione
		movss xmm0, [max]
		mov ecx, [ebp+startt]			;ind2=start
		imul ecx, 4						; solito problema  poi modificare con shift
		add ecx, [edx+dest]				;ind2=start+dest
		push ecx						;  quarto parametro funzione dist

		pushad							;inizio stampa
		push esi
		push dword d1
		call printf
		add esp, 12
		popad

		xor edi, edi
forJ:	cmp edi, [edx+n_centroidi]
		jge endJ						;controllare condizione

		pushad							;inizio stampa
		push edi
		push dword d1
		call printf
		add esp, 12
		popad							;fine stampa

		call distanza
		movss xmm1, [temp]
		comiss xmm1, xmm0
		jge endif
		movss xmm0, xmm1
		mov [eax], edi
endif:	pop ecx
		add ecx, ebx
		push ecx
		jmp forJ
endJ:	pop ecx
		pop ecx
		add ecx, ebx
		push ecx
		mov ecx, [edx+index_columns]
		imul ecx, 4
		add eax, ecx
		jmp forI
endI:	pop	edi		;fine
		pop	esi
		pop	ebx
		mov	esp, ebp
		pop	ebp
		ret


;eax ind
;ebx endj
;ecx temporaneo puntatori
;edx data
;esi i
;edi j

