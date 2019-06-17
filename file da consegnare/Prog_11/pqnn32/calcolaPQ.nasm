extern distanza
extern printf

section .data
max dd 1.79E+38
;f db '%f ', 10, 0
d1 db '%d ', 10, 0
d2 db 'i%d j%d ', 10, 0
;break db 'breakpoint', 10, 0

section .bss
temp: resd 1

section .text

global calcolaPQ

end equ 20
startt equ 16
partition equ 12
data equ 8

source equ 0
dim_source equ 4
index equ 8
dest equ 12
index_rows equ 16
index_columns equ 20
n_centroidi equ 24
d equ 28

; punto1 equ -28
; punto2 equ -24
; dimensione equ -20
; r equ -16

punto1 equ -48
punto2 equ -44
dimensione equ -40
r equ -36

calcolaPQ:
		push ebp
		mov	ebp, esp
		pushad		;inizio
		sub esp, 16
		mov edx, [ebp+data]				;edx=data

		mov ecx, [ebp+partition]		;ind=partition
		sal ecx, 2						; solito problema
		add ecx, [edx+index]			;ind=partition+data->index

		;push temp						;  primo parametro funzione dist
		mov dword [ebp+r], temp


		mov eax, [ebp+end]
		sub eax, [ebp+startt]

		;push eax						;  secondo parametro funzione dist
		mov [ebp+dimensione], eax

		mov eax, [ebp+startt]			;ind1=start
		sal eax, 2						; solito problema  poi modificare con shift
		add eax, [edx+source]			;ind1=start+source

		;push eax						;  terzo parametro funzione dist
		mov [ebp+punto2], eax

		mov ebx, [edx+d]
		sal ebx, 2

		xor esi, esi		;end loop i
forI:	cmp esi, [edx+dim_source]
		jge endI						;controllare condizione
		movss xmm5, [max]
		mov eax, [ebp+startt]			;ind2=start
		sal eax, 2						; solito problema  poi modificare con shift
		add eax, [edx+dest]				;ind2=start+dest

		;push eax						;  quarto parametro funzione dist
		mov [ebp+punto1], eax

		xor edi, edi
forJ:	cmp edi, [edx+n_centroidi]
		jge endJ						;controllare condizione

		; pushad
		; push edi
		; push esi
		; push dword d2
		; call printf
		; add esp, 12
		; popad

		call distanza
		;movss xmm0, [temp]		;istruzione non necessaria, la funzione distanza ha già messo il risultato in xmm0
		
		comiss xmm0, xmm5
		jae endif
		movss xmm5, xmm0

		mov [ecx], edi

endif:	;pop eax
		;add eax, ebx
		;push eax
		add [ebp+punto1], ebx

		inc edi
		jmp forJ
endJ:

		;pop eax
		;add esp, 4	;sostituisce pop eax, forse è meglio
		;pop eax
		;add eax, ebx
		;push eax
		add [ebp+punto2], ebx

		mov eax, [edx+index_columns]
		sal eax, 2
		add ecx, eax

		inc esi
		jmp forI
endI:	add esp, 16
		popad		;fine
		mov	esp, ebp
		pop	ebp
		ret


;ecx ind
;ebx 4*d
;eax temporaneo puntatori e minimo
;edx data
;esi i
;edi j
