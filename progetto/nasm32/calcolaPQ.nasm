section .data
dd max 1.79E+308

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

		mov edi, [ebp+data]				;edi=data
		mov ecx, [ebp+startt]			;ind=start
		imul ecx, [edi+index_columns]	;ind=start*index_columns
		div ecx, [edi+d]				;ind=start*index_columns/d
		add ecx, [edi+index] 			;ind=start*index_columns/d+index
		mov edx, [edi+source]			;ind1=source
		add edx, [ebp+startt]			;ind1=source+start

		xor esi, esi					;i=0
forI:	cmp esi, [edi+dim_source]
		jg endI							;controllare condizione di fine ciclo
		movss xmm2, max
		mov eax, [edi+dest]				;ind2=dest
		add eax, [ebp+startt]			;ind2=dest+start

		xor esp, esp					;j=0
forJ:	cmp esp, [edi+n_centroidi]
		jg endJ							;controllare condizione di fine ciclo
		mov ebx, [edi+end]
		sub ebx, [edi+startt]

		pushad
		push temp
		push ebx
		push eax
		push edx
		call distanza
		;add esp, 16 problema serio
		popad

		; continuare il codice dall'if

		pop	edi		;fine
		pop	esi
		pop	ebx
		mov	esp, ebp
		pop	ebp
		ret