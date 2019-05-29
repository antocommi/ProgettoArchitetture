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
		imul ecx, 4						; solito problema  poi modificare con shift
		add ecx, [edi+index] 			;ind=start*index_columns/d+index
		mov edx, [ebp+startt]			;ind1=start
		imul edx, 4						; solito problema  poi modificare con shift
		add edx, [edi+source]			;ind1=start+source

		xor esi, esi					;i=0
forI:	cmp esi, [edi+dim_source]
		jg endI							;controllare condizione di fine ciclo
		movss xmm2, max
		mov eax, [ebp+startt]			;ind2=start
		imul eax, 4						; solito problema  poi modificare con shift
		add eax, [edi+dest]				;ind2=start+dest

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

		;end if
		
		mov ebx, [edi+d]
		imul ebx, 4
		add eax, ebx
		inc esp
		jmp forJ
endJ:	mov ebx, [edi+index_columns]
		imul ebx, 4
		add ecx, ebx
		mov ebx, [edi+d]
		imul ebx, 4
		add edx, ebx
		inc esi
		jmp forI
endI:	pop	edi		;fine
		pop	esi
		pop	ebx
		mov	esp, ebp
		pop	ebp
		ret

; eax	ind2
; ebx	movimenti
; ecx	ind
; edx	ind1
; esi	i
; edi	data
; ebp
; esp	j